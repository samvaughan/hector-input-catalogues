"""
Take our master region catalogues and remove galaxies which meet various flags.

Do this using RA and DEC limits set in a separate file (resources/RegionInformation/all_regions.csv)

We also remove targets which have bright stars too close to them. The limits are:

* Any star within 5" of the galaxy coordinate
* A 16th magnitude star within the bundle radius
* A 15th magnitude star which is within 2 bundle radii
* A 12th magnitude star or brighter which is within 5 bundle radii
"""
import pandas as pd
import utils
from astroquery.gaia import Gaia
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np

smk = snakemake  # noqa

master_df = pd.read_parquet(smk.input.master_catalogue)

region_information = pd.read_csv(smk.input.region_information, index_col="name")
region_information = region_information.loc[
    region_information.master_region == smk.wildcards.master_region
]

# The largest bundle radius we have- to be conservative
bundle_radius = 20 * u.arcsec

for region, row in region_information.iterrows():
    photometry_source = row.photometry_source

    region_mask = (
        (master_df.RA > row.min_RA)
        & (master_df.RA < row.max_RA)
        & (master_df.DEC > row.min_DEC)
        & (master_df.DEC < row.max_DEC)
    )
    region_catalogue = master_df.loc[region_mask]

    # BAD CLASS 1- too close to stars
    # Now get the stars surrounding our galaxies, with a padding of 10 arcmin
    min_RA = region_catalogue.RA.min() - 10 / 60.0
    max_RA = region_catalogue.RA.max() + 10 / 60.0
    min_DEC = region_catalogue.DEC.min() - 10 / 60.0
    max_DEC = region_catalogue.DEC.max() + 10 / 60.0

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
        dec=region_catalogue.DEC.values * u.degree,
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
    master_df.loc[master_df.ID.isin(bad_class_1_galaxies.ID), "bad_class"] = 1

    # indices_to_remove = bad_class_1_indices
    # # Save our region catalogue
    # postprocessed_region_catalogue = region_catalogue.loc[
    #     ~region_catalogue.index.isin(indices_to_remove)
    # ]
    # postprocessed_region_catalogue.to_csv(output_filename, index=False)
    # print(f"\tSelected {len(postprocessed_region_catalogue)} galaxies for {region}")
