import numpy as np
import astropy.units as u
from scipy.interpolate import interp1d

from .rebin import rebin_1d_trans_box
from .stddata import datadir,dark_sky_Vmag_AB


class Atmosphere:
    '''Atmosphere properties

    This includes the sky transmission and the sky emission.  

    *****

    Sky emission: dark sky emission at CAHA is observational, and obtained from 
        http://www.caha.es/sanchez/sky/
        WARNING: The website claims that the spectra of the sky background are for a 
        PMAS/PPAK fiber (2.7" diam.), normalized to 10-16 Erg s-1 cm-2 A-1. 
        This means the units are 1e-16 erg/s/cm2/A/fiber. However, after a sanity check 
        (comparing CAHA dark sky to Paranal dark sky) it was discovered that the dark
        sky spectrum is already in per arcsec2, not in per fiber. 

    *****

    Sky transmission: The Kelz+2006 paper gives total (telescope+instrument + atmospheric)
    transmission  for a point source at airmass 1.3. To approximately get the atmospheric
    transmission, I divide the total transmission(instrument+telescope + atmosphere) 
    by the telescope+instrument transmission. This gives the approximate atmospheric
    transmission at airmass X = 1.3:
    
    ATM_TRANS(X=1.3) = TRANS_TOTAL / TRANS_TEL_INST
    
    
    since the atmospheric transmission at airmass X approximately goes as 
    ATM_TRANS(X>1) = (ATM_TRANS_ZENITH)^(X>1)
    
    where X is the airmass, I get the atmospheric transmission at zenith as
    ATM_TRANS_ZENITH = (ATM_TRANS(X>1) )^ (1/X)
    
    or in this case:
    
    ATM_TRANS_ZENITH = (ATM_TRANS(X=1.3))^(1/1.3)

    *****

    Seeing calculation is not included.

    The implementation takes one data file from the ESO sky model that
    was taken with standard parameters. For the moon brightness,
    coefficients for moon-target and moon-sun separation are applied.

    The most convenient way to create a CAHA standard atmosphere is
    the function :func:`Atmosphere.caha`.

    Parameters
    ----------

    wavelength : astropy.units.Quantity
        Central wavelength bins for all other arrays [nm]

    zenith_transmission : numpy.ndarray
        Approximate atmospheric transmission for airmass = 1.0 (zenith) 
        absorption lines not included

    dark_sky_flux : astropy.units.Quantity
        Observed dark sky flux in Sanchez et al., PASP 119, 1186  (2007)
        sky spectrum in [photon/(s m² arcsec² nm)]

    sky_mag : astropy.units.Quantity
        Magnitude (per arcsec2) of the sky in AB magnitudes 

    sky_filter_name : str
        For now, must be one of 
        'johnson,b'
        'johnson,v'
        'cousins,r'
        'cousins,i'

    pivot_wave : astropy.units.Quantity
        pivot wavelength of the filter, in which the sky emission was given
        default is the pivot wavelength of the Vband filter (hard-coded) [AA]


    Examples
    --------

    The following example creates a CAHA atmosphere, 
    adjusts the sky emission to match the provided sky magnitude 
    and adds sky transmission. Adding sky transmission is optional:

        >>> atm = Atmosphere.caha()
        >>> sky_emission = atm.emission(20,'johnson,v')
         Multiplying dark sky flux by 1.683
        >>> sky_transmission = atm.transmission(airmass=1.2)
        

    Alternatively, the airmass may be given as zenith angle:

        >>> atm = Atmosphere.caha()
        >>> sky_emission = atm.emission(20,'johnson,v')
         Multiplying dark sky flux by 1.683
        >>> sky_transmission = atm.transmission(airmass=65*u.deg)

    '''

    def __init__(self, wavelength, zenith_transmission,
                 dark_sky_flux,sky_mag,sky_filter_name,pivot_wave):
        self.wavelength = wavelength
        self.zenith_transmission = zenith_transmission
        self.dark_sky_flux = dark_sky_flux
        self.sky_mag = sky_mag
        self.sky_filter_name = sky_filter_name
        self.pivot_wave = pivot_wave


    def transmission(self, airmass=1.0, wavelength=None):
        '''Return the air transmission for a specified airmass.

        Parameters
        ----------

        airmass : float or astropy.units.Quantity
            Airmass to use (default 1.0). The airmass may also be given as
            zenith angle [deg].

        wavelength : astropy.units.Quantity
            Optional wavelength array. If given, the transmission is rebinned
            to this wavelength array, otherwise the unbinned data is returned.

        Returns
        -------
        pysynphot.ArrayBandpass
            ArrayBandpass with transmission coefficients at each wavelength

        '''
        import pysynphot as S
        
        if (isinstance(airmass, u.Quantity)
                and airmass.unit.physical_type == 'angle'):
            airmass = 1.0 / np.cos(airmass)
        res = self.zenith_transmission**airmass
        if wavelength is None:
            return S.ArrayBandpass(self.wavelength,res,name='sky_emission')
        else:
            res = rebin_1d_trans_box(res, self.wavelength, wavelength)
            return S.ArrayBandpass(wavelength,res,name='sky_emission')

    def emission(self,sky_mag,sky_filter_name, wavelength=None):
        '''Return the dark sky emission spectrum, scaled to give integrated surface brightness
        `sky_mag` in filter `sky_filter_name`. 'sky_mag' is in mag/arcsec2.

        the CAHA dark sky spectrum is scaled to the sky_mag in the 
        user selected filter

        Per definition, the zeropoint of an AB mag spectrum is 48.6 in all filters.

        the sky flux F inside of a filter with transmission T (photon-counting transmission curve) is:
        integral(wavelength * T * F) / integral (wavelength * T)

        airmass is not taken into consideration, because there's no moon flux or distance to moon. 


        Parameters
        ----------

        sky_mag : float 
            the magnitude (per arcsec2) of the sky. Dark sky spectrum
            will be scaled to match this magnitude in the provided filter

        sky_filter_name : str 
            the name of the filter in which the sky magnitude is given. 
            must exist in pysynphot. For now, must be one of 
            'johnson,b'
            'johnson,v'
            'cousins,r'
            'cousins,i'

        wavelength : astropy.units.Quantity
            Optional wavelength array. If given, the transmission is rebinned
            to this wavelength array, otherwise the unbinned data is returned.

        Returns
        -------
        pysynphot.ArraySpectrum
            Arrayspectrum with flam as a function of wavelengths
        '''

        import pysynphot as S
        import astropy.constants as const
        
        
        FILTER = S.ObsBandpass(sky_filter_name)
        self.pivot_wave = FILTER.pivot()
        
        dark_flux = self.dark_sky_flux * u.arcsec**2
        
        spec = S.ArraySpectrum(self.wavelength, self.dark_sky_flux,waveunits='angstrom',fluxunits='flam')
        
        obs = S.Observation(spec,FILTER,binset=self.wavelength,force='taper')

        integrated_sky_flux_in_filter = obs.effstim('flam')
        desired_sky_flux_in_filter  = 10**(-0.4*(sky_mag+48.6)) * const.c.to('AA/s').value / FILTER.pivot()**2
        
        scaling_constant = desired_sky_flux_in_filter / integrated_sky_flux_in_filter
        
        print('Multiplying dark sky flux by %.3f'%scaling_constant)
        
        if wavelength is None:
            return S.ArraySpectrum(self.wavelength,self.dark_sky_flux* scaling_constant,waveunits='angstrom',fluxunits='flam') 
        else:
            rebin_flux = rebin_1d_flux_box(spec.flux,spec.wave,wavelength)
            return S.ArraySpectrum(wavelength,rebin_flux,waveunits=wavelength.unit,fluxunits='flam')
        
        
        

    @staticmethod
    def caha():
        '''Standard Calar Alto atmosphere

        This is the atmosphere pieced together from 

        * the Kelz_2006_fig16_ppak_throughput.txt file
        (obtained from Ruben Garcia-Benito via email on November 16, 2020)
        data file 'CAHA_atmo_trans.fits' must exist in  datadir

        * dark sky emission at CAHA 
        (Obtained from http://www.caha.es/sanchez/sky/)
        data file 'dark_sky_emission.fits' must exist in  datadir

        Returns
        -------
        :class:`Atmosphere`
            Atmosphere model created from the standard data files
            included in the package 
        '''

        from astropy.table import QTable
        t = QTable.read(datadir / 'CAHA_atmo_trans.fits')
        s = QTable.read(datadir / 'dark_sky_emission.fits')

        # this will be self.wavelength
        wavelength = t['wavelength'].to('AA').data

        # interpolate dark sky to be on the same wavelength as self.wavelength
        func = interp1d(s['wavelength'].to('AA'),s['dark_sky_flux'].to('erg/(s * cm**2 * AA * arcsec**2)'))
        dark_sky_flux = func(wavelength)
        ###zenith_transmission = func(wavelength)
        
        #sky_filter = Table(data=[[],[]], names=['wavelength','transmission'])
        flux_unit = u.Unit('erg/(s * cm**2 * AA * arcsec**2)')

        
        return Atmosphere(np.copy(wavelength) * u.Angstrom,
                          np.copy(t['zenith_transmission']),
                          np.copy(dark_sky_flux) * flux_unit,
                          dark_sky_Vmag_AB,
                          'johnson,v',
                          5479.35188*u.AA)
#
    
