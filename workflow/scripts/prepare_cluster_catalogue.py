"""
Prepare the Cluster catalogues to be used for Hector target selection

This script does the following things:

* Matches to the redshifts observed with the Hector Redshift Survey
* Adds the stellar mass column, assuming an astropy FlatLambdaCDM cosmology
* Drops galaxies with bad values of stellar mass and/or redshift
"""

import pandas as pd
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
import utils
import pandas_tools as P
import numpy as np


def z_over_sigma_z(x):
    """Take a series of redshifts and find z / sigma(z) so we can remove massive outliers

    Args:
        x (DataFrame): A dataframe with columns 'z' and 'cluster_redshift'

    Returns:
        Series: a series with values z / sigma(z)
    """
    residual = x.z - x.cluster_redshift
    sigma = np.std(residual)
    return residual / sigma


if __name__ == "__main__":

    smk = snakemake
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

    print("Dropping everything without a redshift and which isn't a cluster member")
    good_data_mask = (
        (cluster_catalogue["z"] > 0.0)
        & (cluster_catalogue["logMstarProxy_LS"] > 0.0)
        & ((cluster_catalogue["mem_flag"] == 1))
    )
    cluster_catalogue = cluster_catalogue.loc[good_data_mask]
    print("\tDone")

    # Remove SAMI
    print("Removing SAMI galaxies...")
    sami = P.load_FITS_table_in_pandas(smk.input.SAMI_catalogue)

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
    
    print("Removing a handful of redshift outliers...")
    # Add the cluster redshifts for each object
    cluster_catalogue['cluster_redshift'] = cluster_catalogue.apply(lambda x: cluster_info.loc[x.NearClus, 'z'], axis=1)
    
    # Remove 5 sigma redshift outliers (there's not many)

    redshift_over_sigma_z = cluster_catalogue.groupby('NearClus').apply(z_over_sigma_z).reset_index(level=0, drop=True)
    redshift_over_sigma_z.name = 'z_residual_over_sigma_z'
    cluster_catalogue = cluster_catalogue.join(redshift_over_sigma_z)
    # Now remove everything which has z_over_sigma_z > 5
    z_outlier_mask = cluster_catalogue.loc[:, 'z_residual_over_sigma_z'] < 5
    cluster_catalogue = cluster_catalogue.loc[z_outlier_mask]
    print(f"\tRemoved {np.sum(~z_outlier_mask)} cluster members which are outliers in redshift")
    print("\tDone")

    print("Saving the final catalogue...")
    cluster_catalogue.to_parquet(smk.output.final_cluster_catalogue)
    print("\tDone!")
