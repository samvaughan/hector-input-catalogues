"""
Select standard and guide stars for each of our survey regions in the WAVES regions. For the clusters, we use Matt's star selection.

For WAVES North we use Panstarrs, for WAVES South we use SkyMapper.
"""

import numpy as np
from astroquery.gaia import Gaia
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u
import scipy.spatial


def assign_priorities(df):
    priorities = pd.Series(data=np.zeros(len(df)), name="priority", index=df.index)

    m8 = (df["r_mag"] <= 17) & (df["StandardStar_X_Value"] <= 0.6)
    m7 = (
        (df["r_mag"] > 17) & (df["r_mag"] <= 17.5) & (df["StandardStar_X_Value"] <= 0.6)
    )
    m6 = (
        (df["r_mag"] > 17.5) & (df["r_mag"] <= 18) & (df["StandardStar_X_Value"] <= 0.6)
    )

    m5 = (
        (df["r_mag"] <= 17)
        & (df["StandardStar_X_Value"] > 0.6)
        & (df["StandardStar_X_Value"] <= 1.0)
    )
    m4 = (
        (df["r_mag"] > 17)
        & (df["r_mag"] <= 17.5)
        & (df["StandardStar_X_Value"] > 0.6)
        & (df["StandardStar_X_Value"] <= 1.0)
    )
    m3 = (
        (df["r_mag"] > 17.5)
        & (df["r_mag"] <= 18.0)
        & (df["StandardStar_X_Value"] > 0.6)
        & (df["StandardStar_X_Value"] <= 1.0)
    )

    m2 = (df["r_mag"] <= 17) & (df["StandardStar_X_Value"] > 1)
    m1 = (df["r_mag"] > 17) & (df["StandardStar_X_Value"] > 1)

    priorities.loc[m8] = 8
    priorities.loc[m7] = 7
    priorities.loc[m6] = 6
    priorities.loc[m5] = 5
    priorities.loc[m4] = 4
    priorities.loc[m3] = 3
    priorities.loc[m2] = 2
    priorities.loc[m1] = 1

    return priorities


def standard_star_priority_skymapper(df):
    """Add a priority column for the SkyMapper stars. This priority is minimised for F-type stars (hopefully!)

    Args:
        df (pd.DataFrame): A dataframe containing the magnitudes we need (u, v, g, r, i, z and GAIA BP and RP magnitudes)

    Returns:
        pd.Series: The priority X values
    """

    # These values come from the "Stellar_Colour_Cuts/measure_F_star_colours.py" file
    skymapper_colours = pd.read_csv(smk.input.F_star_colours_SKYMAPPER)
    X = np.sqrt(
        ((df.u_psf - df.v_psf) - skymapper_colours["u_m_v"].values) ** 2
        + ((df.v_psf - df.g_psf) - skymapper_colours["v_m_g"].values) ** 2
        + ((df.g_psf - df.r_psf) - skymapper_colours["g_m_r"].values) ** 2
        + ((df.r_psf - df.i_psf) - skymapper_colours["r_m_i"].values) ** 2
        + ((df.i_psf - df.z_psf) - skymapper_colours["i_m_z"].values) ** 2
        + (
            (df.phot_bp_mean_mag - df.phot_rp_mean_mag)
            - skymapper_colours["Gaia_bp_m_Gaia_rp"].values
        )
        ** 2
    )

    return X


def standard_star_priority_panstarrs(df):
    """Add a priority column for the PANSTARRS stars. This priority is minimised for F-type stars (hopefully!)

    Args:
        df (pd.DataFrame): A dataframe containing the magnitudes we need (g, r, i, z, y and GAIA BP and RP magnitudes)

    Returns:
        pd.Series: The priority X values
    """
    # These values come from the "Stellar_Colour_Cuts/measure_F_star_colours.py" file
    skymapper_colours = pd.read_csv(smk.input.F_star_colours_PANSTARRS)
    X = np.sqrt(
        ((df.z_psf - df.y_psf) - skymapper_colours["z_m_y"].values) ** 2
        + ((df.g_psf - df.r_psf) - skymapper_colours["g_m_r"].values) ** 2
        + ((df.r_psf - df.i_psf) - skymapper_colours["r_m_i"].values) ** 2
        + ((df.i_psf - df.z_psf) - skymapper_colours["i_m_z"].values) ** 2
        + (
            (df.phot_bp_mean_mag - df.phot_rp_mean_mag)
            - skymapper_colours["Gaia_bp_m_Gaia_rp"].values
        )
        ** 2
    )

    return X


def panstarrs_SQL_query(
    min_RA, max_RA, min_DEC, max_DEC, faintest_magnitude, brightest_magnitude, N=500000
):
    panstarrs_query = f"""SELECT TOP {N}
        
        panstarrs.*,
        matchy.source_id,
        matchy.original_ext_source_id,
        main_gaia.ref_epoch,
        main_gaia.ra,
        main_gaia.dec,
        main_gaia.pmra,
        main_gaia.pmdec,
        main_gaia.phot_g_mean_mag,
        main_gaia.phot_bp_mean_mag,
        main_gaia.phot_rp_mean_mag
        
        FROM gaiadr2.panstarrs1_original_valid AS panstarrs

        JOIN gaiadr2.panstarrs1_best_neighbour as matchy
        ON matchy.original_ext_source_id = panstarrs.obj_id
    
        JOIN gaiadr2.gaia_source AS main_gaia ON main_gaia.source_id = matchy.source_id
        
        WHERE (panstarrs.ra > {min_RA})
        AND (panstarrs.ra < {max_RA})
        AND (panstarrs.dec > {min_DEC})
        AND (panstarrs.dec < {max_DEC})
        AND (main_gaia.pmra  > -100)
        AND (main_gaia.pmra < 100)
        AND (main_gaia.pmdec > -100)
        AND (main_gaia.pmdec < 100)
        AND (panstarrs.r_mean_psf_mag > {brightest_magnitude})
        AND (panstarrs.r_mean_psf_mag < {faintest_magnitude})
        AND (panstarrs.g_mean_psf_mag IS NOT NULL)
        AND (panstarrs.r_mean_psf_mag IS NOT NULL)
        AND (panstarrs.i_mean_psf_mag IS NOT NULL)
        AND (panstarrs.z_mean_psf_mag IS NOT NULL)
        AND (panstarrs.y_mean_psf_mag IS NOT NULL)

        ORDER BY r_mean_psf_mag ASC
        """

    return panstarrs_query


def skymapper_SQL_query(
    min_RA, max_RA, min_DEC, max_DEC, faintest_magnitude, brightest_magnitude, N=500000
):
    skymapper_query = f"""SELECT TOP {N}
        
        skymapper.*,
        main_gaia.ref_epoch,
        main_gaia.ra,
        main_gaia.dec,
        main_gaia.pmra,
        main_gaia.pmdec,
        main_gaia.phot_g_mean_mag,
        main_gaia.phot_bp_mean_mag,
        main_gaia.phot_rp_mean_mag
        
        
        FROM external.skymapperdr2_master AS skymapper
        
        JOIN gaiadr3.gaia_source AS main_gaia ON main_gaia.source_id = skymapper.gaia_dr2_id1
        
        WHERE (skymapper.raj2000 > {min_RA})
        AND (skymapper.raj2000 < {max_RA})
        AND (skymapper.dej2000 < {max_DEC})
        AND (skymapper.dej2000 > {min_DEC})
        AND (main_gaia.pmra  > -100)
        AND (main_gaia.pmra < 100)
        AND (main_gaia.pmdec > -100)
        AND (main_gaia.pmdec < 100)
        AND (skymapper.r_psf > {brightest_magnitude})
        AND (skymapper.r_psf < {faintest_magnitude})
        AND (skymapper.flags = 0)
        AND (skymapper.nimaflags = 0)
        AND (skymapper.u_ngood > 0)
        AND (skymapper.g_ngood > 0)
        AND (skymapper.r_ngood > 0)
        AND (skymapper.i_ngood > 0)
        AND (skymapper.z_ngood > 0)

        ORDER BY r_psf ASC
        """

    return skymapper_query


if __name__ == "__main__":
    smk = snakemake  # noqa
    print(f"Region is {smk.wildcards['master_region']}: {smk.wildcards['region_name']}")

    min_RA = smk.params.min_RA - smk.params.pad
    max_RA = smk.params.max_RA + smk.params.pad

    min_DEC = smk.params.min_DEC - smk.params.pad
    max_DEC = smk.params.max_DEC + smk.params.pad

    brightest_magnitude = smk.params.brightest_magnitude
    faintest_magnitude = smk.params.faintest_magnitude

    if smk.wildcards["master_region"] == "WAVES_N":
        query_source = "PANSTARRS"
    elif smk.wildcards["master_region"] == "WAVES_S":
        query_source = "skymapper"
    elif smk.wildcards["master_region"] == "HectorClusters":
        if smk.wildcards["region_name"] in (["A0085", "A0119", "A0151", "A2399"]):
            query_source = "PANSTARRS"
        else:
            query_source = "skymapper"
    else:
        raise NameError("Not sure where to get these stars from!")

    # This is the SQL query we use to select stars from APASS and their match to GAMA.
    # If WAVES North, we use PANSTARRS. Otherwise we use skymapper
    # For the clusters, we use Skymapper for all fields except A0085, A0119 and A2399
    if query_source == "PANSTARRS":
        sql_query = panstarrs_SQL_query(
            min_RA=min_RA,
            max_RA=max_RA,
            min_DEC=min_DEC,
            max_DEC=max_DEC,
            faintest_magnitude=faintest_magnitude,
            brightest_magnitude=brightest_magnitude,
        )
        print("Running the GAIA/PANSTARRS query...")
    elif query_source == "skymapper":
        # Fix an issue we have because the G23 region has a negative RA in our region table
        if smk.wildcards["region_name"] == "G23":
            min_RA += 360
            max_RA += 360
            assert (
                (min_RA > 0) & (min_RA < 360) & (max_RA > 0) & (max_RA < 360)
            ), "The RA constraints must be between 0 and 360"

        sql_query = skymapper_SQL_query(
            min_RA=min_RA,
            max_RA=max_RA,
            min_DEC=min_DEC,
            max_DEC=max_DEC,
            faintest_magnitude=faintest_magnitude,
            brightest_magnitude=brightest_magnitude,
        )
        print("Running the GAIA/SKYMAPPER query...")
    else:
        raise NameError(
            f"The query source must be one of PANSTARRS or skymapper: currently {query_source}"
        )

    # Run the jobs to match to Skymapper/PANSTARRS and GAIA
    job = Gaia.launch_job_async(sql_query)
    r = job.get_results()
    print(f"\tDone! Query selected {len(r)} stars")

    df = r.to_pandas()

    # If we have a panstarrs table, rename one way
    if query_source == "PANSTARRS":
        df.rename(
            columns=dict(
                ra="RA",
                dec="DEC",
                obj_id="ID",
                g_mean_psf_mag="g_psf",
                r_mean_psf_mag="r_psf",
                i_mean_psf_mag="i_psf",
                z_mean_psf_mag="z_psf",
                y_mean_psf_mag="y_psf",
            ),
            inplace=True,
        )
        df.dropna(inplace=True, subset=["g_psf", "r_psf", "i_psf", "z_psf", "y_psf"])
    elif query_source == "skymapper":
        # Otherwise, rename according to the skymapper column names
        df.rename(
            columns=dict(raj2000="RA", dej2000="DEC", object_id="ID"), inplace=True
        )
        df.dropna(inplace=True, subset=["u_psf", "g_psf", "r_psf", "i_psf", "z_psf"])
    else:
        raise NameError(
            f"The query source must be one of PANSTARRS or skymapper: currently {query_source}"
        )

    # Now remove stars which have close companions (i.e. binaries, chance alignments, etc)
    # The max separation is 30 arcseconds, converted to degrees
    max_sep = 30 / 60 / 60
    print(
        f"Removing stars which have companions within {max_sep * 60 *60} arcseconds..."
    )
    X = np.column_stack((df.RA.values, df.DEC.values))
    tree = scipy.spatial.KDTree(X)

    # Query_ball_tree finds all points which are within r of another point
    # It returns a list of lists- for each point, you get the point itself and any neighbours within r
    # So we can then loop through this list and remove stars with more than 1 element in its list
    neighbours = tree.query_ball_tree(tree, r=max_sep)

    stars_to_keep = []
    for neighbour in neighbours:
        if len(neighbour) == 1:
            stars_to_keep.append(neighbour[0])

    N_stars_before_removing_companions = len(df)
    df = df.iloc[stars_to_keep]
    N_stars_after_removing_companions = len(df)
    N_removed_stars = (
        N_stars_before_removing_companions - N_stars_after_removing_companions
    )
    print(f"\tDone. Removed {N_removed_stars} stars")

    print("Accounting for proper motion...")
    # Account for Proper Motion
    # We need a large value of distance in there to avoid an astropy bug- see https://github.com/astropy/astropy/issues/11747
    skycoords = SkyCoord(
        ra=df.RA.values * u.deg,
        dec=df.DEC.values * u.deg,
        pm_ra_cosdec=df.pmra.values * u.mas / u.yr,
        pm_dec=df.pmdec.values * u.mas / u.yr,
        obstime=Time(df.ref_epoch.values, format="jyear"),
        distance=20 * u.kpc,
    )
    # Now update these sky-coords to account for proper motion
    updated_skycoords = skycoords.apply_space_motion(
        Time(smk.params.date_for_observations)
    )

    # Update the RA and DEC values
    df["RA_pmcorr"] = updated_skycoords.ra.value
    df["DEC_pmcorr"] = updated_skycoords.dec.value
    print("\tDone")

    print("Selecting stars and assigning priorities...")
    # Select guide stars
    guide_star_mask = (df.g_psf < 14.5) & (df.g_psf > 14.0)
    guide_stars = df.loc[guide_star_mask].copy()
    hexabundle_stars = df.loc[(~guide_star_mask) & (df.r_psf > 16)].copy()

    # Now select standard stars, sorting by their "priority" values (i.e. colours like an F star)
    if query_source == "PANSTARRS":
        hexabundle_stars["StandardStar_X_Value"] = standard_star_priority_panstarrs(
            hexabundle_stars
        )
    elif query_source == "skymapper":
        hexabundle_stars["StandardStar_X_Value"] = standard_star_priority_skymapper(
            hexabundle_stars
        )
    else:
        raise NameError(
            f"Query_source must be one of PANSTARRS or skymapper; currently {query_source}"
        )

    # Now only keep standard stars with small values of priority
    standard_star_mask = hexabundle_stars["StandardStar_X_Value"] < 1.5
    standard_stars = hexabundle_stars.loc[standard_star_mask].copy()

    # Rename columns to be consistent with other tables.
    column_renamer = dict(
        u_psf="u_mag",
        v_psf="v_mag",
        g_psf="g_mag",
        r_psf="r_mag",
        i_psf="i_mag",
        z_psf="z_mag",
        phot_g_mean_mag="GAIA_g_mag",
        phot_bp_mean_mag="GAIA_bp_mag",
        phot_rp_mean_mag="GAIA_rp_mag",
        pmra="pmRA",
        pmdec="pmDEC",
    )

    # Now set up the column names that we need
    standard_stars.rename(columns=column_renamer, inplace=True)
    guide_stars.rename(columns=column_renamer, inplace=True)

    # Add a couple more
    standard_stars["type"] = 0
    guide_stars["type"] = 2

    # Add the zero columns we expect
    standard_stars.loc[:, ["Mstar", "Re", "z", "GAL_MU_E_R", "y_mag"]] = 0.0

    # And add the priorities
    guide_stars["priority"] = 8
    standard_stars["priority"] = assign_priorities(standard_stars)

    # Finally, subtract 360 from the G23 stars
    if smk.wildcards["region_name"] == "G23":
        guide_stars.RA -= 360
        standard_stars.RA -= 360

    # Only save a subset of the columns
    required_columns_standards = [
        "ID",
        "RA",
        "DEC",
        "g_mag",
        "r_mag",
        "i_mag",
        "z_mag",
        "y_mag",
        "GAIA_g_mag",
        "GAIA_bp_mag",
        "GAIA_rp_mag",
        "Mstar",
        "Re",
        "z",
        "GAL_MU_E_R",
        "pmRA",
        "pmDEC",
        "priority",
        "type",
    ]
    required_columns_guides = ["ID", "RA", "DEC", "r_mag", "type", "pmRA", "pmDEC"]

    standard_stars.loc[:, required_columns_standards].to_csv(
        smk.output.region_standard_star_catalogue, index=False
    )
    guide_stars.loc[:, required_columns_guides].to_csv(
        smk.output.region_guide_star_catalogue, index=False
    )
    print(
        f"\tDone. Selected {len(standard_stars)} standard stars and {len(guide_stars)} guide stars for the region"
    )
