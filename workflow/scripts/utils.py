"""
Useful functions for the catalogue preparation.

Functions are:

* mag_from_flux: calculate magnitudes from a value fo flux (in Janskies)
* apply_Bryant_Mstar_eqtn: Calculates a stellar mass given a redshift, i- and g- band magnitude
* combine_redshift_catalogues: Combines our initial redshift catalogue with the redshifts from the HRS
"""

import numpy as np
import pandas as pd
from numpy.typing import ArrayLike
import typing
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import Table
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib as mpl
import calc_kcor


def mag_from_flux(flux):
    """Make an AB magnitude from a flux in Janskies

    Args:
        flux (float): a flus in Janksies

    Returns:
        float: an AB magnitude
    """

    return -2.5 * np.log10(flux) + 8.90


def apply_Bryant_Mstar_eqtn(redshifts, g_m_i, mag_i, cosmo):
    """Add a stellar mass column using the equation from Bryant et al. 2015

    Args:
        redshifts (array-like): Redshift values
        g_m_i (array-like): g - i colour
        mag_i (array-like): i-band magnitudes in the AB system
        cosmo (astropy.cosmology): An astropy cosmology instance

    Returns:
        array_like: Stellar masses corresponding to the magnitudes in the dataframe, assuming the given cosmology
    """
    # TODO: Fix an error when the redshifts are a pandas series
    dist_mod = cosmo.distmod(redshifts).value

    mstar_effective = (
        -0.4 * mag_i
        + 0.4 * dist_mod
        - np.log10(1.0 + redshifts)
        + (1.2117 - 0.5893 * redshifts)
        + (0.7106 - 0.1467 * redshifts) * g_m_i
    )

    return mstar_effective


def combine_redshift_catalogues(existing_redshifts, observed_redshifts):
    """Combine our catalogue of existing redshifts with the redshifts from the Hector Redshift Survey observations. Only keep existing redshifts which have 'usez == 1'

    Args:
        existing_redshifts (pd.DataFrame): Existing redshift dataframe
        observed_redshifts (pd.DataFrame): Hector Redshift Survey dataframe

    Returns:
        _type_: _description_
    """

    # Only keep things which have 'usez=1' for the existing redshift catalogue and rename some columns
    existing_redshifts = existing_redshifts.loc[existing_redshifts.usez == 1]
    existing_redshifts = existing_redshifts.loc[:, ["ra", "dec", "z_helio"]].rename(
        dict(ra="RA", dec="DEC", z_helio="Z"), axis=1
    )
    observed_redshifts = observed_redshifts.loc[:, ["RA", "DEC", "Z"]]

    # Combine the two together
    all_redshifts = pd.concat((existing_redshifts, observed_redshifts)).reset_index()

    return all_redshifts


def match_catalogues(
    catalogue_A: tuple[ArrayLike, ArrayLike], catalogue_B: tuple[ArrayLike, ArrayLike]
) -> typing.Tuple[ArrayLike, ArrayLike, ArrayLike]:
    """Match two catalogues in redshift using Astropy and return the idx, distance and 3D distance lists. Note! The longest catalogue should be given first- the output distance vector is the same length as catalogue B.

    Args:
        catalogue_A (tuple[ArrayLike, ArrayLike]): A two component tuple of RA and Dec vectors for catalogue A
        catalogue_B (tuple[ArrayLike, ArrayLike]): A two component tuple of RA and Dec vectors for catalogue B

    Returns:
        typing.Tuple[ArrayLike, ArrayLike, ArrayLike]: A tuple of the idx, 2D (on-sky) distance and 3D distance between each match.
    """

    catalogue_A_RA, catalogue_A_DEC = catalogue_A
    catalogue_B_RA, catalogue_B_DEC = catalogue_B

    cat_A = SkyCoord(ra=catalogue_A_RA * u.degree, dec=catalogue_A_DEC * u.degree)
    cat_B = SkyCoord(
        ra=catalogue_B_RA * u.degree,
        dec=catalogue_B_DEC * u.degree,
    )
    idx, distance, threed_distance = cat_B.match_to_catalog_sky(cat_A)

    return idx, distance, threed_distance


def select_MW_analogues(df):
    mwa_mask = (
        (10.5 < df["Mstar"])
        & (df["Mstar"] < 10.9)
        & (0.875 < df["g_m_i"])
        & (df["g_m_i"] < 1.275)
    )

    return mwa_mask


def select_one_to_two_r200(df):
    mask = (df["r_on_rtwo"] > 1) & (df["r_on_rtwo"] < 2)

    return mask


def select_within_r200(df):
    mask = df["r_on_rtwo"] < 1

    return mask


def select_outside_2r200(df):
    mask = df["r_on_rtwo"] > 2

    return mask


def select_edge_on_wind_galaxies(df):
    wind_mask = (1 - df["Ellipticity_r"] < 0.5) & (df["g_m_i"] < 0.6)

    return wind_mask


def add_MWA_edge_on_priority(df):
    """Select Galaxies from Jesse's MWA selection or Brent's edge-on wind galaxy selection

    Args:
        df (pd.DataFrame): a DataFrame with columns 'Mstar', 'g_m_i', 'Ellipticity_r' and 'r_mag'

    Returns:
        pd.Series: A boolean mask which is True for galaxies which pass the selection
    """

    mwa_mask = select_MW_analogues(df)
    edge_on_wind_mask = select_edge_on_wind_galaxies(df)

    overall_mask = mwa_mask | edge_on_wind_mask

    return overall_mask


def low_mass_galaxy_mask(df, mass_name="Mstar", g_m_i_name="g_m_i"):
    """Select galaxies with a low stellar mass and red colours.

    Args:
        df (pd.DataFrame): a DataFrame with columns 'Mstar' and 'g_m_i'

    Returns:
        pd.Series: A boolean mask which is True for galaxies which pass the selection
    """

    low_mass_mask = (df[mass_name] < 9.5) & (df[g_m_i_name] > 0.5)

    return low_mass_mask


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


def k_corrected_r_band_abs_mag(df, cosmo):
    # Calculate the absolute R-band magnitude, including a K correction
    k_cors = calc_kcor.calc_kcor("r", df.z, "g - r", df.mag_g - df.mag_r)
    abs_mag_r = (
        df.mag_r.values * u.mag - cosmo.distmod(df.z.values) - k_cors.values * u.mag
    ).value

    return abs_mag_r


def get_membership_prob(samples, data):
    membership_prob = 1 - np.exp(
        samples.loc[:, "log_Pr[1]":f"log_Pr[{data['N']}]"].mean()
    )
    return membership_prob


def plot_red_sequence(samples, data, mean_abs_mag):
    """Plot the colour-magnitude diagram of the clusters and show the straight line fits and the membership probabilities

    Args:
        samples (_type_): _description_
        data (_type_): _description_
        mean_abs_mag (_type_): _description_

    Returns:
        fig, ax:
    """

    # Make some straight lines
    xx = np.linspace(-4, 4)
    # # Blue Cloud
    # yy_blue = (
    #     xx[:, None] * np.array(samples["slopes[1]"])[None, :]
    #     + np.array(samples["intercepts[1]"])[None, :]
    # )
    # red sequence
    yy_red = (
        xx[:, None] * np.array(samples["slopes[2]"])[None, :]
        + np.array(samples["intercepts[2]"])[None, :]
    )

    # Make the probabilities
    Red_Sequence_member_probability = get_membership_prob(samples, data)

    # Make some plots
    # Firstly, the red sequence membership probability on the color/magnitude diagram
    # plt.style.use("publication")
    fig, ax = plt.subplots(figsize=(14, 6))
    ax.scatter(
        data["x"] + mean_abs_mag,
        data["y"],
        c=Red_Sequence_member_probability.values,
        alpha=0.2,
        cmap="RdYlBu_r",
        rasterized=True,
        s=100,
    )

    # Make the colorbar, without any alpha
    cmap = mpl.cm.RdYlBu_r
    norm = mpl.colors.Normalize(vmin=0, vmax=1)
    fig.colorbar(
        mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
        ax=ax,
        label=r"$p(\mathrm{red\,\,sequence\,\,member})$",
    )

    # Add on some straight line fits
    # ax.plot(xx + mean_abs_mag, yy_blue.mean(axis=1), c='k', alpha=1, linestyle='dotted', linewidth=1.5)
    ax.plot(
        xx + mean_abs_mag,
        yy_red.mean(axis=1),
        c="k",
        alpha=1,
        linestyle="solid",
        linewidth=3.5,
    )
    ax.plot(
        xx + mean_abs_mag,
        yy_red.mean(axis=1) - 2 * samples["scatter[2]"].mean(),
        c="0.2",
        linestyle="dashed",
        linewidth=1.5,
    )
    ax.plot(
        xx + mean_abs_mag,
        yy_red.mean(axis=1) + 2 * samples["scatter[2]"].mean(),
        c="0.2",
        linestyle="dashed",
        linewidth=1.5,
    )

    ax.set_xlabel("Absolute $r$-band magnitude (mag)")
    ax.set_ylabel("$g - r$ colour (mag)")
    ax.set_ylim(0.0, 1.2)

    return fig, ax


def load_FITS_table_in_pandas(filename):
    hdu = fits.open(filename)
    t = Table(hdu[1].data)
    return t.to_pandas()


def skymapper_SQL_query(
    min_RA,
    max_RA,
    min_DEC,
    max_DEC,
    faintest_magnitude,
    brightest_magnitude,
    N=500000,
    require_all_bands=True,
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
        AND (skymapper.r_psf > {brightest_magnitude})
        AND (skymapper.r_psf < {faintest_magnitude})
        """
    if require_all_bands:
        skymapper_query += """AND (main_gaia.pmra  > -100)
        AND (main_gaia.pmra < 100)
        AND (main_gaia.pmdec > -100)
        AND (main_gaia.pmdec < 100)
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


def panstarrs_SQL_query(
    min_RA,
    max_RA,
    min_DEC,
    max_DEC,
    faintest_magnitude,
    brightest_magnitude,
    N=500000,
    require_all_bands=True,
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
        AND (panstarrs.r_mean_psf_mag > {brightest_magnitude})
        AND (panstarrs.r_mean_psf_mag < {faintest_magnitude})
        """
    if require_all_bands:
        panstarrs_query += """AND (main_gaia.pmra  > -100)
        AND (main_gaia.pmra < 100)
        AND (main_gaia.pmdec > -100)
        AND (main_gaia.pmdec < 100)
        AND (panstarrs.g_mean_psf_mag IS NOT NULL)
        AND (panstarrs.r_mean_psf_mag IS NOT NULL)
        AND (panstarrs.i_mean_psf_mag IS NOT NULL)
        AND (panstarrs.z_mean_psf_mag IS NOT NULL)
        AND (panstarrs.y_mean_psf_mag IS NOT NULL)

        ORDER BY r_mean_psf_mag ASC
        """

    return panstarrs_query
