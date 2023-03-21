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

if __name__ == "__main__":
    
    smk = snakemake
    sep_constraint_arcsec = smk.params.sep_constraint_arcsec
    
    # Make the desired cosmology
    cosmology = FlatLambdaCDM(H0=70, Om0=0.3)

    cluster_catalogue = pd.read_parquet(smk.input.input_cluster_catalogue)
    
    # Get the redshift catalogues
    print("Combining the redshift catalogues...")
    existing_redshifts = pd.read_parquet(smk.input.exisiting_redshift_catalogue).reset_index()
    observed_redshifts = pd.read_parquet(smk.input.observed_redshift_catalogue)
    all_redshifts = utils.combine_redshift_catalogues(existing_redshifts, observed_redshifts)
    print("\tDone")
    
    # Perform the matching
    print("Matching redshifts to the cluster catalogues...")
    redshift_cat = (all_redshifts.RA.values, all_redshifts.DEC.values)
    photom_cat = (cluster_catalogue.RA.values, cluster_catalogue.dec.values)
    idx, distance, threed_distance = utils.match_catalogues(redshift_cat, photom_cat)
    # Only keep things within sep arcseconds of one another
    match_mask = distance < sep_constraint_arcsec * u.arcsec
    
    cluster_catalogue.loc[match_mask, 'z'] = all_redshifts.loc[
        all_redshifts.index[idx[match_mask]], "Z"
    ].values
    print("\tDone")
    
    cluster_catalogue["logMstarProxy_LS"] = utils.apply_Bryant_Mstar_eqtn(
        cluster_catalogue.z.values, cluster_catalogue.LS_FAKE_SDSS_gmi.values, cluster_catalogue.LS_FAKE_SDSS_i.values, cosmology
    )
    
    print("Dropping everything without a redshift...")
    good_data_mask = (cluster_catalogue['z'] > 0.0) & (cluster_catalogue['logMstarProxy_LS'] > 0.0)
    cluster_catalogue = cluster_catalogue.loc[good_data_mask]
    print("\tDone")
    
    print("Saving the final catalogue...")
    cluster_catalogue.to_parquet(smk.output.final_cluster_catalogue)
    print("\tDone!")
    
