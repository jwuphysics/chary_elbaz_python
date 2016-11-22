# John F. Wu
# 2016-05-16
"""
A simple port of the Chary & Elbaz (2001) SED fitting code to python. 

Note that the numerical integration routine we use (scipy.integrate.simps) 
is different from the IDL function (INT_TABULATED), which uses a five-point
Newton-Cotes method (as opposed to the three-point Newton-Cotes formula in
our python implementation). I found a ~0.03% difference between luminosity 
values in my own use cases.
"""

import os
import numpy as np
import astropy.units as u
from astropy.cosmology import FlatLambdaCDM
from scipy.io import readsav
from scipy.integrate import simps

CHARY_ELBAZ_DIR = '/home/john/Research/alma_cycle2_analysis/code/misc/chary_elbaz/'


def quick_fit_SED(z, S_nu, wavelength, use_IRAS=True):
    """
    Essentially a wrapper for the fit_SED function above. S_nu and 
    wavelength are input as unitless scalars but assumed to be in units 
    of uJy and um, respectively. Assumes a concordance cosmology.
    """

    return fit_SED(z, S_nu * u.uJy, wavelength * u.um, use_IRAS=use_IRAS)


@u.quantity_input(S_nu=u.uJy, wavelength=u.um)
def fit_SED(z, S_nu, wavelength, H0=70., Om0=0.3, resolution=3., use_IRAS=True):
    """
    Given a redshift of source and measured flux density at a specified 
    wavelength, returns the total infrared luminosity by fitting to a 
    library of template SEDs. See Chary & Elbaz (2001) for details. This
    program is modeled after the IDL program chary_elbaz.pro.


    Parameters
    ----------
    z (float)          : redshift
    S_nu (float)       : flux density at the specified wavelength in uJy
    wavelength (float) : filter wavelength of measured flux density in
                         microns (um)

    H0 (float)         : Hubble constant; defaults to 70 km s^-1 Mpc^-1
    Om0 (float)        : matter density; defaults to 0.3
    resolution (float) : instrument spectral resolution; defaults to 3.
    Sanders_LIR (bool) : return Sanders L_IR (see Chary & Elbaz 2001 
                         eq. 1) or otherwise the true L_IR ; defaults 
                         to True.

    Returns
    -------
    L_IR : (float) : total infrared luminosity, in units of solar luminosity
    """

    S_nu = S_nu.value
    wavelength = wavelength.value

    # Read in Chary & Elbaz (2001) save file
    ce = readsav(CHARY_ELBAZ_DIR + 'chary_elbaz.save')

    # establish cosmology
    cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)
    D_L = cosmo.luminosity_distance(z).value   # Mpc

    # read in wavelengths and redshift them
    # shape : (1366,)
    wavelengths = ce['lambda']
    wavelengths_observed = wavelengths * (1 + z)

    # read in templates, redshift, and convert to flux density
    # shape : (1366, 105)
    SEDs = ce['nuLnuinLsun'] 
    SEDs_observed = np.array([(1 + z) * f / (1e-32 * 4 * np.pi * 
        (D_L * 3.0856e22)**2 * (3e14 / wavelengths) / 3.826e26) \
        for f in SEDs.T]).T
    
    # find edges of instrument filter based on resolution
    wavemin = wavelength * (1 - 1. / (2 * resolution))
    wavemax = wavelength * (1 + 1. / (2 * resolution))


    # convert templates into what the instrument sees
    # shape : (105,)
    SEDs_instrument = np.mean(SEDs_observed[(wavelengths_observed > wavemin) &
                                            (wavelengths_observed < wavemax)],
                              axis = 0)

    # index of best match
    ind = np.argmin(np.abs(S_nu - SEDs_instrument))

    # return Chary & Elbaz L_IR (Sanders) if use_IRAS is False
    if use_IRAS == False: 
        L_IR = ce['LIR_Sanders'][ind]
    else:
        # shape : (1366,)
        flux_densities = SEDs_observed[:, ind] * S_nu / SEDs_instrument[ind]

        # used for calculating L_IR (Sanders & Mirabel 1996)
        f_12  = iras(12,  wavelengths, flux_densities / (1 + z))
        f_25  = iras(25,  wavelengths, flux_densities / (1 + z))
        f_60  = iras(60,  wavelengths, flux_densities / (1 + z))
        f_100 = iras(100, wavelengths, flux_densities / (1 + z))

        # Sanders L_IR via IRAS fluxes (Chary & Elbaz 2001, eq 1)
        L_IR = 1.8e-14 * 1e-6 * (13.48 * f_12 + 5.16 * f_25 + 2.58 * f_60 + 
            f_100) * 4 * np.pi * (D_L * 3.0856e22)**2 / 3.826e26

    return L_IR

def iras(iras_filter, wavelengths, flux_densities):
    """
    Fulfills the duty of the iras12.pro, iras25.pro, etc. IDL programs.

    Parameters
    ----------
    iras_filter (int)      : the wavelength of the IRAS filter; must be 12, 
                             25, 60, or 100 (microns)
    wavelengths (array)    : the wavelengths in microns (1366,)

    flux_densities (array) : the flux densities in uJy (1366,)

    Returns
    -------
    f_iras (float) : the integrated flux coming through the IRAS filter 


    """

    if iras_filter not in [12, 25, 60, 100]:
        import sys
        sys.exit('Enter a valid IRAS filter.')

    # shape : (21,)
    iras_wave, iras_resp = \
        np.genfromtxt(CHARY_ELBAZ_DIR + 'iras{}.txt'.format(str(iras_filter)),
                      skip_header=3, unpack=True)

    iras_flux_densities = np.interp(iras_wave, wavelengths, flux_densities)

    fct_numerator = iras_flux_densities * iras_resp / iras_wave**2
    int_numerator = simps(fct_numerator, iras_wave)

    fct_denominator = iras_resp / (iras_wave * iras_filter)
    int_denominator = simps(fct_denominator, iras_wave)

    f_iras = int_numerator / int_denominator

    return f_iras
