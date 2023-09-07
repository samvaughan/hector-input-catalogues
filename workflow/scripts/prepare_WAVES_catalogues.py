"""
Prepare the WAVES input catalogues we've been given to be used for Hector target selection.

This script does the following things:
 
* Removes any objects classed as "star"
* Matches the photmetric catalogues to our redshifts (both existing and observed with the HRS)
* Adds the correct magnitude columns
* Adds the stellar mass column, assuming an astropy FlatLambdaCDM cosmology
* Removes SAMI galaxies from WAVES North
* Adds in the bad_class column and removes stars too close to galaxies
"""

import pandas as pd
import astropy.units as u
import numpy as np
import utils
from astropy.cosmology import FlatLambdaCDM
from astroquery.gaia import Gaia
from astropy.coordinates import SkyCoord

if __name__ == "__main__":
    smk = snakemake  # noqa
    sep_constraint_arcsec = smk.params.sep_constraint_arcsec

    # Make the desired cosmology
    cosmology = FlatLambdaCDM(H0=70, Om0=0.3)

    # Read in the existing redshifts and our new observations
    print("Combining the redshifts...")
    existing_redshifts = pd.read_parquet(
        smk.input.exisiting_redshift_catalogue
    ).reset_index()
    observed_redshifts = pd.read_parquet(smk.input.observed_redshift_catalogue)
    all_redshifts = utils.combine_redshift_catalogues(
        existing_redshifts, observed_redshifts
    )
    print("\tDone")

    # Now load the photometry from WAVES
    print("Combining the WAVES photometry...")
    waves_catalogue_S = pd.read_parquet(smk.input.WAVES_S_input)
    waves_catalogue_N = pd.read_parquet(smk.input.WAVES_N_input)
    all_photometry = pd.concat((waves_catalogue_S, waves_catalogue_N))
    print("\tDone")

    # Remove things classed as "star"
    print("Removing stars...")
    to_remove = all_photometry["class"] == "star"
    all_photometry = all_photometry.loc[~to_remove]
    print("\tDone")

    # Fold the stuff with a negative RA to be around ~330-360 again
    # We need to do this because the redshift catalogue we're matching to runs between 0 and 360 not -180 and 180.
    all_photometry.loc[all_photometry.RAmax < 0, "RAmax"] = (
        all_photometry.loc[all_photometry.RAmax < 0, "RAmax"] + 360.0
    )

    # Do the matching
    print("Matching redshifts to the WAVES catalogues...")
    redshift_cat = (all_redshifts.RA.values, all_redshifts.DEC.values)
    photom_cat = (all_photometry.RAmax.values, all_photometry.Decmax.values)
    idx, distance, threed_distance = utils.match_catalogues(redshift_cat, photom_cat)

    # Only keep things within sep arcseconds of one another
    match_mask = distance < sep_constraint_arcsec * u.arcsec
    print("\tDone")

    # Now drop things without a redshift
    print("Dropping everything without a redshift...")
    # Use the match_mask as a new column and add the redshifts
    all_photometry["has_redshift"] = match_mask
    all_photometry["z"] = np.nan
    all_photometry.loc[match_mask, "z"] = all_redshifts.loc[
        all_redshifts.index[idx[match_mask]], "Z"
    ].values
    all_photometry = all_photometry.loc[match_mask]
    print("\tDone")

    print("Making the magnitude columns...")
    # Now make the magnitude columns
    flux_columns = [
        "flux_ut",
        "flux_gt",
        "flux_rt",
        "flux_it",
        "flux_Zt",
        "flux_Yt",
        "flux_Jt",
        "flux_Ht",
        "flux_Kt",
    ]
    new_mag_columns = [
        "mag_ut",
        "mag_gt",
        "mag_rt",
        "mag_it",
        "mag_Zt",
        "mag_Yt",
        "mag_Jt",
        "mag_Ht",
        "mag_Kt",
    ]
    # Add the magnitude columns
    # Drop any NaNs
    all_photometry = all_photometry.dropna(subset=flux_columns)
    for new_col, flux_col in zip(new_mag_columns, flux_columns):
        all_photometry[new_col] = utils.mag_from_flux(all_photometry[flux_col])

    # Drop the mag_app_ columns since they're pretty useless
    cols_to_drop = [
        "mag_app_ut",
        "mag_app_gt",
        "mag_app_rt",
        "mag_app_it",
        "mag_app_Zt",
        "mag_app_Yt",
        "mag_app_Jt",
        "mag_app_Ht",
        "mag_app_Kt",
    ]
    all_photometry = all_photometry.drop(cols_to_drop, axis=1)
    print("\tDone")

    print("Adding the stellar masses...")
    all_photometry["Mstar"] = utils.apply_Bryant_Mstar_eqtn(
        all_photometry.z.values,
        (all_photometry.mag_gt - all_photometry.mag_it).values,
        all_photometry.mag_it.values,
        cosmology,
    )
    print("\tDone")

    print("Adding the ellipticity...")
    # Make an ellipticity column from the axis ratio
    all_photometry["ellipticity"] = 1 - all_photometry["axrat"]
    print("\tDone")

    # Move things back to have negative RAs
    all_photometry.loc[all_photometry.RAmax > 300, "RAmax"] -= 360

    # Add a bad_class column
    all_photometry["bad_class"] = 0
    """
    A galaxy with bad_class = 1 has:
    * Any star within 5" of the galaxy coordinate
    * A 16th magnitude star within the bundle radius
    * A 15th magnitude star which is within 2 bundle radii
    * A 12th magnitude star or brighter which is within 5 bundle radii
    """

    region_information = pd.read_csv(smk.input.region_information, index_col="name")
    region_information = region_information.loc[
        region_information.master_region.isin(["WAVES_N", "WAVES_S"])
    ]

    # The largest bundle radius we have- to be conservative
    bundle_radius = 20 * u.arcsec
    all_photometry["bad_class"] = 0

    for region, row in region_information.iterrows():
        photometry_source = row.photometry_source

        region_mask = (
            (all_photometry.RAmax > row.min_RA)
            & (all_photometry.RAmax < row.max_RA)
            & (all_photometry.Decmax > row.min_DEC)
            & (all_photometry.Decmax < row.max_DEC)
        )
        assert (
            region_mask.sum() > 0
        ), "Looks like our region mask doesn't match any galaxies!"
        region_catalogue = all_photometry.loc[region_mask]

        # BAD CLASS 1- too close to stars
        # Now get the stars surrounding our galaxies, with a padding of 10 arcmin
        min_RA = region_catalogue.RAmax.min() - 10 / 60.0
        max_RA = region_catalogue.RAmax.max() + 10 / 60.0
        min_DEC = region_catalogue.Decmax.min() - 10 / 60.0
        max_DEC = region_catalogue.Decmax.max() + 10 / 60.0

        if region == "G23":
            min_RA += 360
            max_RA += 360

        if photometry_source == "PANSTARRS":
            sql_query_function = utils.panstarrs_SQL_query
            r_mag_name = "r_mean_psf_mag"
        elif photometry_source == "skymapper":
            sql_query_function = utils.skymapper_SQL_query
            r_mag_name = "r_psf"
        else:
            raise NameError(f"Unknown photometry source! {photometry_source}")

        sql_query = sql_query_function(
            min_RA,
            max_RA,
            min_DEC,
            max_DEC,
            faintest_magnitude=18,
            brightest_magnitude=-5,
            require_all_bands=False,
        )
        print(f"{region}: Getting nearby stars...")
        job = Gaia.launch_job_async(sql_query)
        r = job.get_results()
        print(f"\tDone! Query selected {len(r)} stars")

        stars = r.to_pandas()

        # Now use Astropy to match the two catalogues
        galaxy_coords = SkyCoord(
            ra=region_catalogue.RAmax.values * u.degree,
            dec=region_catalogue.Decmax.values * u.degree,
        )
        star_coords = SkyCoord(
            ra=stars.ra.values * u.degree, dec=stars.dec.values * u.degree
        )
        idx, d2d, d3d = galaxy_coords.match_to_catalog_sky(star_coords)

        # Now work through our classifications
        bad_class_1_indices = []
        print("\tChecking for galaxies near bright stars...")

        # First, any star which is within 2" of the galaxy centre
        mask_1 = d2d < 2 * u.arcsec
        # print(
        #     f"\t\tRemoving {mask_1.sum()} galaxies which have a star within 2`` of their centres"
        # )
        bad_class_1_indices.extend(region_catalogue.index[mask_1])

        # Now any star brighter than 16th mag within the bundle radius
        star_coords = SkyCoord(
            ra=stars.loc[stars[r_mag_name] < 16, "ra"].values * u.degree,
            dec=stars.loc[stars[r_mag_name] < 16, "dec"].values * u.degree,
        )
        idx, d2d, d3d = galaxy_coords.match_to_catalog_sky(star_coords)
        mask_2 = d2d < bundle_radius
        # print(
        #     f"\t\tRemoving {mask_2.sum()} galaxies with stars brighter than 16th mag within the bundle radius"
        # )
        bad_class_1_indices.extend(region_catalogue.index[mask_2])

        # Now any star brighter than 15th mag within 2 bundle radiii
        star_coords = SkyCoord(
            ra=stars.loc[stars[r_mag_name] < 15, "ra"].values * u.degree,
            dec=stars.loc[stars[r_mag_name] < 15, "dec"].values * u.degree,
        )
        idx, d2d, d3d = galaxy_coords.match_to_catalog_sky(star_coords)
        mask_3 = d2d < 2 * bundle_radius
        # print(
        #     f"\t\tRemoving {mask_3.sum()} galaixes with stars brighter than 15th mag within two bundle radii"
        # )
        bad_class_1_indices.extend(region_catalogue.index[mask_3])

        # Now any star brighter than 12th mag within 5 bundle radiii
        star_coords = SkyCoord(
            ra=stars.loc[stars[r_mag_name] < 12, "ra"].values * u.degree,
            dec=stars.loc[stars[r_mag_name] < 12, "dec"].values * u.degree,
        )
        if len(star_coords) > 0:
            idx, d2d, d3d = galaxy_coords.match_to_catalog_sky(star_coords)
            mask_4 = d2d < 5 * bundle_radius
            # print(
            #     f"\t\tRemoving {mask_4.sum()} galaxies with stars brighter than 12th mag within five bundle radii"
            # )
            bad_class_1_indices.extend(region_catalogue.index[mask_3])
        else:
            pass

        bad_class_1_indices = np.unique(bad_class_1_indices)
        print(f"\tRemoved {len(bad_class_1_indices)} galaxies as bad_class 1")

        # Get our removed galaxies and flip the 'bad_class' in the master catalogue
        bad_class_1_galaxies = region_catalogue.loc[
            region_catalogue.index.isin(bad_class_1_indices)
        ]
        all_photometry.loc[
            all_photometry.uberID.isin(bad_class_1_galaxies.uberID), "bad_class"
        ] = 1

    # Split the catalogues up again
    all_photometry["uberID"] = all_photometry.apply(lambda x: f"W{x.uberID}", axis=1)
    final_WAVES_S_catalogue = all_photometry.loc[all_photometry.Decmax < -15]
    waves_N_catalogue_with_SAMI = all_photometry.loc[all_photometry.Decmax > -15]

    print("Removing SAMI galaxies from WAVES North...")
    # Now remove the SAMI galaxies from WAVES_N
    sami = utils.load_FITS_table_in_pandas(smk.input.SAMI_catalogue)

    sami_catalogue = (sami.RA.values, sami.DEC.values)
    waves_N_catalogue = (
        waves_N_catalogue_with_SAMI.RAmax.values,
        waves_N_catalogue_with_SAMI.Decmax.values,
    )
    idx, d2d, d3d = utils.match_catalogues(sami_catalogue, waves_N_catalogue)

    max_sep = sep_constraint_arcsec * u.arcsec
    sep_constraint = d2d < max_sep

    final_WAVES_N_catalogue = waves_N_catalogue_with_SAMI.loc[~sep_constraint]

    print(
        f"\tRemoved {np.sum(sep_constraint)} galaxies from the WAVES North catalogue which are also in the SAMI catalogue"
    )
    print("\tDone")

    # Add in the approximate surface brightness values
    final_WAVES_N_catalogue["approximate_SB_r"] = final_WAVES_N_catalogue[
        "mag_rt"
    ] + 2.5 * np.log10(
        np.pi * final_WAVES_N_catalogue["R50"] * final_WAVES_N_catalogue["R50"]
    )
    final_WAVES_S_catalogue["approximate_SB_r"] = final_WAVES_S_catalogue[
        "mag_rt"
    ] + 2.5 * np.log10(
        np.pi * final_WAVES_S_catalogue["R50"] * final_WAVES_S_catalogue["R50"]
    )

    print("Saving the final catalogues...")
    final_WAVES_N_catalogue.to_parquet(smk.output.final_WAVES_N_catalogue)
    final_WAVES_S_catalogue.to_parquet(smk.output.final_WAVES_S_catalogue)
    print("\tDone!")
