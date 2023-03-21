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


def mag_from_flux(flux):
    """Make an AB magnitude from a flux in Janskies

    Args:
        flux (float): a flus in Janksies

    Returns:
        float: an AB magnitude
    """

    return -2.5 * np.log10(flux) + 8.90


def apply_Bryant_Mstar_eqtn(redshifts, mag_g, mag_i, cosmo):
    """Add a stellar mass column using the equation from Bryant et al. 2015

    Args:
        redshifts (array-like): Redshift values
        mag_g (array-like): g-band magnitudes in the AB system
        mag_i (array-like): i-band magnitudes in the AB system
        cosmo (astropy.cosmology): An astropy cosmology instance

    Returns:
        array_like: Stellar masses corresponding to the magnitudes in the dataframe, assuming the given cosmology
    """

    dist_mod = cosmo.distmod(redshifts).value

    mstar_effective = (
        -0.4 * mag_i
        + 0.4 * dist_mod
        - np.log10(1.0 + redshifts)
        + (1.2117 - 0.5893 * redshifts)
        + (0.7106 - 0.1467 * redshifts) * (mag_g - mag_i)
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


def match_catalogues(catalogue_A: tuple[ArrayLike, ArrayLike], catalogue_B: tuple[ArrayLike, ArrayLike]) -> typing.Tuple[ArrayLike, ArrayLike, ArrayLike]:
    """Match two catalogues in redshift using Astropy and return the idx, distance and 3D distance lists.

    Args:
        catalogue_A (tuple[ArrayLike, ArrayLike]): A two component tuple of RA and Dec vectors for catalogue A
        catalogue_B (tuple[ArrayLike, ArrayLike]): A two component tuple of RA and Dec vectors for catalogue B

    Returns:
        typing.Tuple[ArrayLike, ArrayLike, ArrayLike]: A tuple of the idx, 2D (on-sky) distance and 3D distance between each match.
    """
    
    catalogue_A_RA, catalogue_A_DEC = catalogue_A
    catalogue_B_RA, catalogue_B_DEC = catalogue_B
    
    cat_A = SkyCoord(
        ra=catalogue_A_RA * u.degree, dec=catalogue_A_DEC * u.degree
    )
    cat_B = SkyCoord(
        ra=catalogue_B_RA * u.degree,
        dec=catalogue_B_DEC * u.degree,
    )
    idx, distance, threed_distance = cat_B.match_to_catalog_sky(cat_A)

    return idx, distance, threed_distance