"""
Prepare the WAVES input catalogues we've been given to be used for Hector target selection.

This script does the following things:
 
* Removes any objects classed as "star"
* Matches the photmetric catalogues to our redshifts (both existing and observed with the HRS)
* Adds the correct magnitude columns
* Adds the stellar mass column, assuming an astropy FlatLambdaCDM cosmology
* Removes SAMI galaxies from WAVES North
"""

import pandas as pd
import astropy.units as u
import numpy as np
import pandas_tools as P
from astropy.cosmology import FlatLambdaCDM
import utils


if __name__ == "__main__":

    smk = snakemake
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
    all_photometry['ellipticity'] = 1 - all_photometry['axrat']
    print("\tDone")

    # Split the catalogues up again
    final_WAVES_S_catalogue = all_photometry.loc[all_photometry.Decmax < -15]
    waves_N_catalogue_with_SAMI = all_photometry.loc[all_photometry.Decmax > -15]

    print("Removing SAMI galaxies from WAVES North...")
    # Now remove the SAMI galaxies from WAVES_N
    sami = P.load_FITS_table_in_pandas(smk.input.SAMI_catalogue)

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

    # Finally, move the things we've put around 350 degrees back to be negative in RA
    final_WAVES_S_catalogue.loc[final_WAVES_S_catalogue.RAmax > 300, "RAmax"] -= 360

    print("Saving the final catalogues...")
    final_WAVES_N_catalogue.to_parquet(smk.output.final_WAVES_N_catalogue)
    final_WAVES_S_catalogue.to_parquet(smk.output.final_WAVES_S_catalogue)
    print("\tDone!")
