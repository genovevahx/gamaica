# Idea and some text taken from
# https://github.com/NASA-Planetary-Science/sbpy/blob/master/sbpy/units/core.py
# Copyright (c) 2017-2019, sbpy team (3-Clause BSD License)
"""Vega unit core module

This package defines Vega-based magnitude systems, and a units
equivalency function for converting to/from flux density.

To use these units in the top-level `astropy.units` namespace, you
need to import the module::

    >>> import astropy.units as u
    >>> import etc.units
    >>> u.Unit('mag(VEGA)')
    Unit("mag(VEGA)")

"""

import numpy as np
import astropy.units as u

VEGA = u.def_unit(['VEGA', 'VEGAflux'],
                  doc='Spectral flux density of Vega.')

VEGAmag = u.MagUnit(VEGA)
VEGAmag.__doc__ = "Vega-based magnitude: Vega is 0 mag at all wavelengths"

u.add_enabled_units(VEGA)


def spectral_density(wfb):
    """Flux density equivalencies, extended for Vega-based magnitude systems.

    This uses the default Vega spectrum that comes for 4MOST in the
    Instrument Control Repository (ICR).

    Vega is assumed to have an apparent magnitude of 0 in the
    ``VEGAmag`` system.

    See https://github.com/astropy/astropy/issues/10821 for a
    discussion about using `u.spectral_density` instead.

    Parameters
    ----------
    wfb : `astropy.units.Quantity`
        Wavelength, of the corresponding flux density being converted.


    Returns
    -------
    equiv : list
        List of equivalencies.


    Examples
    --------
    Monochromatic flux density:

        >>> import astropy.units as u
        >>> from etc.units import spectral_density, VEGAmag
        >>> m = 10 * VEGAmag
        >>> fluxd = m.to('erg/(s cm² nm)', spectral_density(5557.5*u.AA))
        >>> print(f'{fluxd:.3e}')
        3.461e-12 erg / (cm2 nm s)
        >>> fluxd = m.to('ph/(s cm² nm)', spectral_density(5557.5*u.AA))
        >>> print(f'{fluxd:.4f}')
        0.9682 ph / (cm2 nm s)

    """
    global _vega
    if _vega is None:
        _vega = VegaSpectrum()

    flux0 = _vega(wfb)
    return [
        (flux0.unit, VEGA,
         lambda f_phys: f_phys / flux0,
         lambda f_vega: f_vega * flux0),
        (u.Unit('ph/(s cm² nm)'), VEGA,
         lambda f_phys: f_phys / flux0.to('ph/(s cm² nm)',
                                          u.spectral_density(wfb)),
         lambda f_phys: f_phys * flux0.to('ph/(s cm² nm)',
                                          u.spectral_density(wfb)))
    ] + u.spectral_density(wfb)


class VegaSpectrum:
    def __init__(self, name='alpha_lyr'):
        from .spectrum import _read_spectrum_from_fits
        self.wavelength, self.flux = _read_spectrum_from_fits(name)

    def __call__(self, wavelength, unit=None):
        return np.interp(wavelength, self.wavelength, self.flux)


_vega = None
