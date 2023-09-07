"""
Prepare the Cluster catalogues to be used for Hector target selection **including foreground galaxies along the line of sight**

This script does the following things:

* Matches to the redshifts observed with the Hector Redshift Survey
* Adds the stellar mass column, assuming an astropy FlatLambdaCDM cosmology
* Drops galaxies with bad values of stellar mass and/or redshift
"""

import pandas as pd
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
import utils
import numpy as np
from cmdstanpy import CmdStanModel
from astroquery.gaia import Gaia
from astropy.coordinates import SkyCoord


if __name__ == "__main__":
    smk = snakemake  # noqa
    sep_constraint_arcsec = smk.params.sep_constraint_arcsec

    # Make the desired cosmology
    cosmology = FlatLambdaCDM(H0=70, Om0=0.3)

    cluster_catalogue = pd.read_parquet(smk.input.input_cluster_catalogue)
    cluster_info = pd.read_table(
        smk.input.cluster_information, delim_whitespace=True, index_col=0
    )

    # Spoke to Matt on 23/03/23: He will deal with the redshift matching, and I'll stop doing it here.

    # # Get the redshift catalogues
    # print("Combining the redshift catalogues...")
    # existing_redshifts = pd.read_parquet(
    #     smk.input.exisiting_redshift_catalogue
    # ).reset_index()
    # observed_redshifts = pd.read_parquet(smk.input.observed_redshift_catalogue)
    # all_redshifts = utils.combine_redshift_catalogues(
    #     existing_redshifts, observed_redshifts
    # )
    # print("\tDone")

    # # Perform the matching
    # print("Matching redshifts to the cluster catalogues...")
    # redshift_cat = (all_redshifts.RA.values, all_redshifts.DEC.values)
    # photom_cat = (cluster_catalogue.RA.values, cluster_catalogue.dec.values)
    # idx, distance, threed_distance = utils.match_catalogues(redshift_cat, photom_cat)
    # # Only keep things within sep arcseconds of one another
    # match_mask = distance < sep_constraint_arcsec * u.arcsec

    # cluster_catalogue.loc[match_mask, "z"] = all_redshifts.loc[
    #     all_redshifts.index[idx[match_mask]], "Z"
    # ].values
    # print("\tDone")

    print("Making the Stellar Mass and Elliptpcity columns...")
    cluster_catalogue["logMstarProxy_LS"] = utils.apply_Bryant_Mstar_eqtn(
        cluster_catalogue.z.values,
        cluster_catalogue.LS_FAKE_SDSS_gmi.values,
        cluster_catalogue.LS_FAKE_SDSS_i.values,
        cosmology,
    )
    cluster_catalogue["ellipticity"] = 1 - cluster_catalogue["B_on_A"]
    print("\tDone")

    # print("Removing a handful of redshift outliers...")
    # Add the cluster redshifts for each object
    cluster_catalogue["cluster_redshift"] = cluster_catalogue.apply(
        lambda x: cluster_info.loc[x.NearClus, "z"], axis=1
    )

    print(
        "Dropping everything without a redshift. Keeping Cluster members as well as galaxies which are within 2.5 r200 of the BCG, are in front of the cluster, and fall into our low-mass three-repeats selection"
    )
    good_data_mask_cluster_members = (
        (cluster_catalogue["z"] > 0.0)
        & (cluster_catalogue["logMstarProxy_LS"] > 0.0)
        & ((cluster_catalogue["mem_flag"] == 1))
        & ((cluster_catalogue["r_on_rtwo"] <= 2.5))
    )
    good_data_mask_foreground = (
        (cluster_catalogue["z"] > 0.0)
        & (cluster_catalogue["logMstarProxy_LS"] > 0.0)
        & (cluster_catalogue["mem_flag"] == 0)
        & (cluster_catalogue["r_on_rtwo"] < 2.5)
        & (cluster_catalogue["logMstarProxy_LS"] < 9.5)
        & (cluster_catalogue["z"] < cluster_catalogue["cluster_redshift"])
    )

    good_data_mask = good_data_mask_foreground | good_data_mask_cluster_members
    cluster_catalogue = cluster_catalogue.loc[good_data_mask]
    print("\tDone")

    # Fix the cluster redshift column for non members
    # Otherwise the foreground things almost never pass the target selection
    cluster_catalogue.loc[
        good_data_mask_foreground, "cluster_redshift"
    ] = cluster_catalogue.loc[good_data_mask_foreground, "z"]

    # Remove SAMI
    print("Removing SAMI galaxies...")
    sami = utils.load_FITS_table_in_pandas(smk.input.SAMI_catalogue)

    sami_catalogue = (sami.RA.values, sami.DEC.values)
    cluster_catalogue_RA_DEC = (
        cluster_catalogue.RA.values,
        cluster_catalogue.dec.values,
    )
    idx, d2d, d3d = utils.match_catalogues(sami_catalogue, cluster_catalogue_RA_DEC)

    max_sep = sep_constraint_arcsec * u.arcsec
    sep_constraint = d2d < max_sep

    cluster_catalogue = cluster_catalogue.loc[~sep_constraint]
    print(f"\tRemoved {np.sum(sep_constraint)} SAMI galaxies")
    print("\tDone")

    # Fit the red sequence
    print("Fitting the red sequence...")
    # first get absolute magnitudes
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    cluster_catalogue["Absolute_mag_R"] = utils.k_corrected_r_band_abs_mag(
        cluster_catalogue, cosmo
    )

    # Centre these values around the mean and make the y values
    mean_abs_mag = np.mean(
        cluster_catalogue.loc[good_data_mask_cluster_members, "Absolute_mag_R"]
    )
    x_all = (
        cluster_catalogue.loc[good_data_mask_cluster_members, "Absolute_mag_R"]
        - mean_abs_mag
    )
    y_all = (
        cluster_catalogue.loc[good_data_mask_cluster_members, "mag_g"]
        - cluster_catalogue.loc[good_data_mask_cluster_members, "mag_r"]
    )

    # Get rid of some galaxies which are very red
    x = x_all.loc[y_all < 1.05]
    y = y_all.loc[y_all < 1.05]
    N = len(x)
    data = dict(N=N, x=x, y=y)

    # Fit the stan model
    sm = CmdStanModel(stan_file=smk.input.stan_file)
    fit = sm.sample(data=data)

    # Get the results as a pandas dataframe
    samples = fit.draws_pd()

    fig, ax = utils.plot_red_sequence(samples, data, mean_abs_mag)
    fig.savefig(smk.output.cluster_RS_plot, bbox_inches="tight")

    # Get the membership probabilities and save the table
    Red_Sequence_member_probability = utils.get_membership_prob(samples, data)
    Red_Sequence_member_probability.index = data["x"].index

    # Using the fact that the index is now correct
    cluster_catalogue["RS_member_probability"] = 0
    cluster_catalogue["RS_member_probability"] = Red_Sequence_member_probability

    # Tweak the fact that some very red things have a smaller probability of being a red sequence member than they should
    cluster_catalogue.loc[
        (good_data_mask_cluster_members)
        & (cluster_catalogue.mag_g - cluster_catalogue.mag_r > 1.0),
        "RS_member_probability",
    ] = 1

    cluster_catalogue["RS_member_from_mm"] = (
        cluster_catalogue["RS_member_probability"] > 0.4
    )

    # Add the column where we take everything above the red sequence minus 2 sigma of scatter to be a red sequence member
    cluster_catalogue["RS_member_2sig_scatter"] = 0.0
    cluster_catalogue.loc[
        good_data_mask_cluster_members, "RS_member_2sig_scatter"
    ] = data["y"] > (
        data["x"] * samples["slopes[2]"].mean()
        + samples["intercepts[2]"].mean()
        - 2 * samples["scatter[2]"].mean()
    )
    cluster_catalogue["RS_member_2sig_scatter"] = cluster_catalogue[
        "RS_member_2sig_scatter"
    ].astype(bool)

    # Add the approximate surface brightnesses
    cluster_catalogue["approximate_SB_r"] = cluster_catalogue["mag_r"] + 2.5 * (
        np.log10(np.pi) + 2 * np.log10(cluster_catalogue["Re"])
    )

    # Add the bad_class attributes
    """
    A galaxy with bad_class = 1 has:
    * Any star within 5" of the galaxy coordinate
    * A 16th magnitude star within the bundle radius
    * A 15th magnitude star which is within 2 bundle radii
    * A 12th magnitude star or brighter which is within 5 bundle radii
    """

    region_information = pd.read_csv(smk.input.region_information, index_col="name")
    region_information = region_information.loc[
        region_information.master_region == "HectorClusters"
    ]

    # The largest bundle radius we have- to be conservative
    bundle_radius = 20 * u.arcsec
    cluster_catalogue["bad_class"] = 0

    # Match with the bad_class attributes from our table
    # We only have visual classifications for ~3k galaxies
    # So assume that everything is good unless it is flagged as bad
    badclass_df = pd.read_csv(smk.input.badclass_votes)
    merged = pd.merge(
        cluster_catalogue,
        badclass_df.loc[:, ["ID", "BAD_CLASS"]],
        left_on="ls_id",
        right_on="ID",
        how="left",
    )
    # merged now has our badclass categories for everything we've looked at
    # using the HTS app, and NaNs for everything we haven't
    # Replace these NaNs with the 0s from the cluster catalogue
    # And then set the 'bad_class' attribute equal to this
    cluster_catalogue["bad_class"] = (
        merged["BAD_CLASS"].fillna(merged["bad_class"]).astype(int)
    ).values
    # Now go through and do the 'near star' classes using the criteria
    # from SAMI
    for region, row in region_information.iterrows():
        photometry_source = row.photometry_source

        region_mask = (
            (cluster_catalogue.RA > row.min_RA)
            & (cluster_catalogue.RA < row.max_RA)
            & (cluster_catalogue.dec > row.min_DEC)
            & (cluster_catalogue.dec < row.max_DEC)
        )
        region_catalogue = cluster_catalogue.loc[region_mask]

        # BAD CLASS 1- too close to stars
        # Now get the stars surrounding our galaxies, with a padding of 10 arcmin
        min_RA = region_catalogue.RA.min() - 10 / 60.0
        max_RA = region_catalogue.RA.max() + 10 / 60.0
        min_DEC = region_catalogue.dec.min() - 10 / 60.0
        max_DEC = region_catalogue.dec.max() + 10 / 60.0

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
            ra=region_catalogue.RA.values * u.degree,
            dec=region_catalogue.dec.values * u.degree,
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

        cluster_catalogue.loc[
            cluster_catalogue.ls_id.isin(bad_class_1_galaxies.ls_id), "bad_class"
        ] = 1

    print("\tDone!")
    cluster_catalogue["ls_id"] = cluster_catalogue.apply(
        lambda x: f"C{x.ls_id}", axis=1
    )
    print("Saving the final catalogue...")
    cluster_catalogue.to_parquet(smk.output.final_cluster_catalogue)
    print("\tDone!")
