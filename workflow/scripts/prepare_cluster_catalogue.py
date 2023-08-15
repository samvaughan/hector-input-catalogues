"""
Prepare the Cluster catalogues to be used for Hector target selection

This script does the following things:

* Matches to the redshifts observed with the Hector Redshift Survey
* Adds the stellar mass column, assuming an astropy FlatLambdaCDM cosmology
* Drops galaxies with bad values of stellar mass and/or redshift
* Fits the red sequence to get red sequence membership probabilities/flags
"""

import pandas as pd
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
import utils
import pandas_tools as P
import numpy as np
from cmdstanpy import CmdStanModel


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
        & ((cluster_catalogue["r_on_rtwo"] <= 2.5))
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

    # Add the cluster redshifts for each object
    cluster_catalogue["cluster_redshift"] = cluster_catalogue.apply(
        lambda x: cluster_info.loc[x.NearClus, "z"], axis=1
    )

    # Fit the red sequence
    print("Fitting the red sequence...")
    # first get absolute magnitudes
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    cluster_catalogue["Absolute_mag_R"] = utils.k_corrected_r_band_abs_mag(
        cluster_catalogue, cosmo
    )

    # Centre these values around the mean and make the y values
    mean_abs_mag = np.mean(cluster_catalogue.Absolute_mag_R)
    x_all = cluster_catalogue.Absolute_mag_R - mean_abs_mag
    y_all = cluster_catalogue.mag_g - cluster_catalogue.mag_r

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
    cluster_catalogue["RS_member_probability"] = Red_Sequence_member_probability

    # Tweak the fact that some very red things have a smaller probability of being a red sequence member than they should
    cluster_catalogue.loc[
        cluster_catalogue.mag_g - cluster_catalogue.mag_r > 1.0, "RS_member_probability"
    ] = 1

    cluster_catalogue["RS_member_from_mm"] = (
        cluster_catalogue["RS_member_probability"] > 0.4
    )

    # Add the column where we take everything above the red sequence minus 2 sigma of scatter to be a red sequence member
    cluster_catalogue["RS_member_2sig_scatter"] = data["y"] > (
        data["x"] * samples["slopes[2]"].mean()
        + samples["intercepts[2]"].mean()
        - 2 * samples["scatter[2]"].mean()
    )
    cluster_catalogue["RS_member_2sig_scatter"] = cluster_catalogue[
        "RS_member_2sig_scatter"
    ].astype(bool)
    
    cluster_catalogue['approximate_SB_r'] = cluster_catalogue['mag_r'] + 2.5 * (np.log10(np.pi) + 2 * np.log10(cluster_catalogue['Re']))

    print("\tDone!")

    print("Saving the final catalogue...")
    cluster_catalogue.to_parquet(smk.output.final_cluster_catalogue)
    print("\tDone!")
