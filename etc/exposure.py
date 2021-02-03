import numpy as np
import astropy.units as u
from astropy.table import QTable

from .rebin import rebin_1d_flux_box,rebin_1d_trans_box

class Exposure:
    '''Simulate exposures

    After initialization and filling with object and background (sky)
    flux, the basic properties for a specified exposure time can be
    calculated.

    Parameters
    ----------
    spectrograph : :class:`.Gamaica`
        Spectrograph object, with attributes wavelength, readout_noise,
        dark_flux,, and others.

    Examples
    --------

    Do a full ETC run

        >>> from etc.exposure import Exposure
        >>> from etc.instrument import Gamaica
        >>> from etc.atmosphere import Atmosphere
        >>> from etc.seds import SEDtemplate
        >>> from etc.functions import save_and_see
        >>> 
        >>> 
        >>>## instrument
        >>>instrument = Gamaica(1200)
        >>>
        >>>## atmosphere
        >>>atm = Atmosphere.caha()
        >>>sky_emission = atm.emission(21.5,'johnson,v')
        >>>sky_transmission = atm.transmission(airmass=1.2)
        >>>
        >>># target to observe
        >>>HII_region = SEDtemplate('orion')
        >>>HII_region.scale_to_mag(20.0,'johnson,v')
        >>>
        >>>
        >>># set up exposure
        >>>exp = Exposure(instrument)
        >>>
        >>># add object to observe
        >>>exp.add_object(HII_region)
        >>>
        >>># add sky emission, optional
        >>>exp.add_sky_emission(sky_emission)
        >>>
        >>># add atmospheric transmission (optional)
        >>>exp.add_sky_transmission(sky_transmission)
        >>>
        >>># expose
        >>>exptime = 60. # seconds
        >>>result = exp.expose(exptime)
        >>>
        >>># save result, look at plots
        >>>save_and_see(result)
    '''

    def __init__(self, spectrograph):
        self.spectrograph = spectrograph
        #self.object_flux = np.zeros(self.spectrograph.wavelength.shape) * (
        #    u.electron * u.s**-1)
        #self.sky_flux = np.zeros(self.spectrograph.wavelength.shape) * (
        #    u.electron * u.s**-1)

        self.nobj_A = np.zeros(self.spectrograph.wavelength_perpix.shape) * ( u.photon/u.second/u.m**2/u.AA)
            
        self.nobj_fib = np.zeros(self.spectrograph.wavelength_perpix.shape) *( u.photon/u.second/u.AA)
        self.nobj_pix = np.zeros(self.spectrograph.wavelength_perpix.shape) *( u.photon/u.second)
        self.nobj_photlam = np.zeros(self.spectrograph.wavelength_perpix.shape) *( u.photon/u.second/u.cm**2/u.AA)
        self.nobj_fiber_frac = 1.0
        
        self.nsky_A = np.zeros(self.spectrograph.wavelength_perpix.shape) * ( u.photon/u.second/u.m**2/u.AA)
        self.nsky_fib = np.zeros(self.spectrograph.wavelength_perpix.shape) *( u.photon/u.second/u.AA)
        self.nsky_pix = np.zeros(self.spectrograph.wavelength_perpix.shape) *( u.photon/u.second)
        self.nsky_photlam = np.zeros(self.spectrograph.wavelength_perpix.shape) *(  u.photon/u.second/u.cm**2/u.AA)
        
        
    def add_object(self, sed, wavelength=None):
        '''Add some object flux

        Parameters
        ----------
        sed : :class:`SEDtemplate`
            object flux density in [flam] as a function of wavelength in [A]

        wavelength : :class:`astropy.units.Quantity`
            Wavelength bins for the flux. If not given, then the 
            flam is not rebinned and the wavelengths of each bin must 
            correspond to the wavelength bins of the GAMAICA transmission
        '''

        
        ## conversions
        sed.sed.convert('photlam')
        ''' this conversion seems to strip sed.sed.flux from its units
        but not sed.sed.wave! Inconvenient.
        Probably best to replace pysynphot with own code.
        '''


        # rebin incoming sed to instrument dispersion
        # rebin doesn't handle units, resulting in errors.
        # avoid by stripping arrays of units
        if isinstance(sed.sed.flux,np.ndarray):
            flux_perpix = rebin_1d_flux_box(sed.sed.flux.value,sed.sed.wave.value,self.wavelength.value)* u.photon/u.second/u.cm**2/u.AA
        elif isinstance(sed.sed.flux,u.quantity.Quantity):
            flux_perpix = rebin_1d_flux_box(sed.sed.flux.value,sed.sed.wave.value,self.wavelength.value)* u.photon/u.second/u.cm**2/u.AA
        else:
            print('Cannot add object. ****not implemented****')
            return
        nobj_A = flux_perpix.to('ph s**-1 m**-2 AA**-1')
        nobj_fib = nobj_A * self.spectrograph.A_tel * self.spectrograph.fib_area.value


        ## this formula below for nobj_pix needs checking. I assume that 1 AA = 1/ldisp pixels. Correct?
        nobj_pix = nobj_fib * self.spectrograph.ldisp
        
        ## find index of wavelength element closest to pivot position
        idx = abs(self.wavelength.to('AA').value - sed.pivot_wave).argmin()
        
        print(f'\nAt pivot wavelength {self.wavelength[idx]:.0f} of filter:')
        print(f'Object flam: {sed.sed.flux[idx]:.5e} erg / (Angstrom cm2 s)' )
        print(f'Object photons/m2: {nobj_A[idx]:.5e} ')
        print(f'Spatial sampling: {self.spectrograph.fib_area:.3f}')
        #print(f'fraction of star/fiber: {(self.spectrograph.fib_diam/fwhm)**2)}')
        print(f'Object photons/fiber:  {nobj_fib[idx]:.5e}')
        print(f'Object photons/pixel:  {nobj_pix[idx]:.5e}')

        self.nobj_photlam += flux_perpix
        self.nobj_A += nobj_A
        self.nobj_fib += nobj_fib
        #self.fiber_frac = (fib/fwhm)**2
        self.nobj_pix += nobj_pix
        
        if wavelength is not None:
            print('rebinning of target flux ***Not implemented yet***')
        #    flux = self.spectrograph.rebin(wavelength, flux)
        

    def add_sky_emission(self, sed, wavelength=None):
        '''Add some sky flux

        Parameters
        ----------
        sed : :class:`pysynphot.spectrum.CompositeSourceSpectrum` or
                     `pysynphot.spectrum.SourceSpectrum`
            flux density of sky in erg/s/cm2/AA.
            This is actually sky surface brightness. The "per arcsec2" 
            is implied but not explicit, otherwise the unit conversions
            don't work.

        wavelength : :class:`astropy.units.Quantity`
            Wavelength bins for the flux. If not given, then the 
            flam is not rebinned and the wavelengths of each bin must 
            correspond to the wavelength bins of the GAMAICA transmission

        '''
        
        sed.convert('photlam')
        ''' this conversion seems to strip sed.sed.flux from its units
        but not sed.sed.wave! Inconvenient.'''

        
        # rebin incoming sky sed to instrument dispersion
        # rebin doesn't handle units, resulting in errors.
        # avoid by stripping arrays of units
        if isinstance(sed.flux,np.ndarray):
            flux_perpix = rebin_1d_flux_box(sed.flux.value,sed.wave.value,self.wavelength.value)* u.photon/u.second/u.cm**2/u.AA
        elif isinstance(sed.flux,u.quantity.Quantity):
            flux_perpix = rebin_1d_flux_box(sed.flux.value,sed.wave.value,self.wavelength.value)* u.photon/u.second/u.cm**2/u.AA
        else:
            print('Cannot add sky emission. ****not implemented****')
            return
        
        nsky_A = flux_perpix.to('ph s**-1 m**-2 AA**-1')
        nsky_fib = nsky_A * self.spectrograph.A_tel * self.spectrograph.fib_area.value
        

        ## this formula below for nobj_pix needs checking. I assume that 1 AA = 1/ldisp pixels. Correct?
        nsky_pix = nsky_fib * self.spectrograph.ldisp
        
        ## print result for default Vband for now
        ## find index of wavelength element closest to pivot position
        idx = abs(self.wavelength.to('AA').value - 5479.35188).argmin()

        print(f'\nAt pivot wavelength 5479 Angstrom of Vband:')
        print(f'Sky flam: {sed.flux[idx]:.5e} erg / (Angstrom cm2 s)' )
        print(f'Sky photons/m2: {nsky_A[idx]:.5e} ')
        print(f'Spatial sampling: {self.spectrograph.fib_area:.3f}')
        #print(f'fraction of sky/fiber: {(self.spectrograph.fib_diam/fwhm)**2)}')
        print(f'Sky photons/fiber:  {nsky_fib[idx]:.5e}')
        print(f'Sky photons/pixel:  {nsky_pix[idx]:.5e}')


        self.nsky_photlam += flux_perpix
        self.nsky_A += nsky_A
        self.nsky_fib += nsky_fib
        self.nsky_pix += nsky_pix
        
        
        if wavelength is not None:
            print('rebinning of sky flux ***Not implemented yet***')
        #    flux = self.spectrograph.rebin(wavelength, flux)
        

    def add_sky_transmission(self, trans, wavelength=None):
        '''Add some sky flux

        Parameters
        ----------
        trans : :class:`pysynphot.spectrum.CompositeSpectralElement` or
                       `pysynphot.spectrum.ArraySpectralElement`
            atmospheric transmission at the user-specified airmass
            

        wavelength : :class:`astropy.units.Quantity`
            Wavelength bins for the transmission. If not given, then it
            is not rebinned and the wavelengths of each bin must 
            correspond to the wavelength bins of the GAMAICA wavelength vector

        '''
        
        ## TODO:
        # 1. decide on rebinning
        # 2. get values at pivot wavelength
        
        

        # rebin incoming sky sed to instrument dispersion
        trans_perpix = rebin_1d_trans_box(trans.throughput,trans.wave.value,self.wavelength.value)
        
        if wavelength is not None:
            print('rebinning of sky flux ***Not implemented yet***')
        #    flux = self.spectrograph.rebin(wavelength, flux)
        
        # add atmospheric transmission to tel+instr transmission
        self.spectrograph.efficiency_perpix *=  trans_perpix


    def expose(self, texp, nexp=1, xbin=4,ybin=4):
        '''Calculate the result of one exposure

        WARNING: Seeing calculation is not included. This means all of the
        target's flux is assumed to enter the fiber.


        This basically scales the object and sky fluxes by the
        exposure time, adds dark and readout noise, and calculates the
        signal-to-noise ratio. The result is returned as a
        :class:`astropy.table.Qtable`.

        The noise is calculated by the formula
        :math:`\\sqrt{(N_{obj} + N_{sky} + N_{dark}) ⋅ t_{exp} \
              + N_{ron}^2 ⋅ n_{exp}}`,
        where all values are in electrons.


        Parameters
        ----------
        texp : :class:`astropy.units.Quantity`
            Total exposure time [s]

        nexp : int
            Number of exposures. Defaults to 1.

        xbin : int
            Binning in x direction

        ybin : int
            Binning in y direction

        Returns
        -------
        :class:`astropy.table.QTable` (if no columns were specified)
            Resulting table with the following columns:

              * ``wavelength``: Wavelength [AA]
              * ``ldisp``: Dispersion per pixel [AA]
              * ``exptime``: exposure time [s]
              * ``object``: Object count [e-]
              * ``sky``: Sky background count [e-]
              * ``dark``: CCD dark current [e-]
              * ``ron``: CCD readout noise [e-]
              * ``noise``: Noise count [e-]

        '''

        texp *= u.second
        
        #--- compute number of photons per pixel ------------------------
        # counts per pixel over total exposure: projected fiber image on (un)binned
        # CCD, system efficiency, reciprocal linear dispersion, exposure time

        ## binning. Currently no binning
        nobj_pix = self.nobj_pix * xbin *ybin /16.0

        ## tel+instr+atmo throughput (Q.E. included)
        nobj_pix = nobj_pix * self.spectrograph.efficiency_perpix * texp
        
        ## binning. Currently no binning
        nsky_pix = self.nsky_pix * xbin *ybin /16.0
        
        ## tel+instr+atmo throughput (Q.E. included)
        nsky_pix = nsky_pix * self.spectrograph.efficiency_perpix * texp

        # summarize noise contributions

        #shotnoise contributions per pixel
        shotn_obj  = np.sqrt(nobj_pix)
        shotn_sky  = np.sqrt(nsky_pix)

        #dark current
        ndark_pix   = self.spectrograph.dark_current * xbin * ybin * texp
        shotn_dark   = np.sqrt(ndark_pix)
        
        # noise total
        noise_pix = np.sqrt((shotn_obj**2).to_value(u.electron) +
                            (shotn_sky**2).to_value(u.electron) +
                            (shotn_dark**2).to_value(u.electron) +
                            self.spectrograph.ron.to_value(u.electron)**2 ) *u.electron



        #--- Signal-to-Noise per wavelength bin ------------
        nwbin = 4/ybin	# number of pixels to sum-up for one wavelength bin
        
        signal = nwbin * nobj_pix
        backgr = nwbin * nsky_pix
        #cal_noise =  backgr  * cal_limit
        
        
        #noise  = sqrt( nwbin * (noise_pix)**2 + (cal_noise)**2 )
        noise  = np.sqrt( nexp * nwbin * (noise_pix)**2 )
        s2n    = (nexp*signal)/noise


        return QTable({'wavelength': self.wavelength,
                       'ldisp':np.repeat(self.spectrograph.ldisp,
                                         len(self.wavelength)),
                       #'gain': np.repeat(self.spectrograph.gain,
                       #                  len(self.wavelength)),
                       'exptime':np.repeat(texp,
                                         len(self.wavelength)),
                       'object': nobj_pix,
                       'sky': nsky_pix,
                       'dark': np.repeat(ndark_pix,len(self.wavelength)),
                       'ron': np.repeat(np.sqrt(nexp) * self.spectrograph.ron, len(self.wavelength)),
                       'noise': noise,
                       'snr': s2n * u.pixel**-1})

    @property
    def wavelength(self):
        return self.spectrograph.wavelength_perpix
