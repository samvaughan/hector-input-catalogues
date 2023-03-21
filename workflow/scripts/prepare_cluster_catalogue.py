"""
Prepare the Cluster catalogues to be used for Hector target selection

This script does the following things:

* Matches to the redshifts observed with the Hector Redshift Survey
* Adds the stellar mass column, assuming an astropy FlatLambdaCDM cosmology
* Drops galaxies with bad values of stellar mass and/or redshift
"""

import pandas as pd
import numpy as np
import astropy.units as u
import pandas_tools as P
from astropy.cosmology import FlatLambdaCDM
import utils 

if __name__ == "__main__":
    
    smk = snakemake
    sep_constraint_arcsec = smk.params.sep_constraint_arcsec
    
    cluster_catalogue = pd.read_parquet(smk.input.input_cluster_catalogue)
    
    # Get the redshift catalogues
    print("Combining the redshift catalogues...")
    existing_redshifts = pd.read_parquet(smk.input.exisiting_redshift_catalogue).reset_index()
    observed_redshifts = pd.read_parquet(smk.input.observed_redshift_catalogue)
    all_redshifts = utils.combine_redshift_catalogues(existing_redshifts, observed_redshifts)
    print("\tDone")
    
    # Perform the matching
    redshift_cat = (all_redshifts.RA.values, all_redshifts.DEC.values)
    photom_cat = (cluster_catalogue.RA.values, cluster_catalogue.dec.values)
    idx, distance, threed_distance = photom_cat.match_to_catalog_sky(redshift_cat, photom_cat)
    # Only keep things within sep arcseconds of one another
    match_mask = distance < sep_constraint_arcsec * u.arcsec
    print("\tDone")
