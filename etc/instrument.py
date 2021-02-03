import numpy as np
import astropy.units as u
import astropy.table
from astropy.table import QTable


from .rebin import rebin_1d_trans_box
import etc.stddata as defaults

       
class Gamaica:

    '''Combines Telescope and Instrument properties

    This class describes the light path through the common parts of
    the telescope, i.e. from the mirror, through the fiber, to the detector.


    Parameters
    ----------

    A_tel : astropy.units.Quantity
        Effective area of the telescope [m²]

    fib_area : astropy.unit.Quantity
        Area of the fiber, projected to the sky. [arcsec2]

    fib_diam : astropy.unit.Quantity
        Diameter of fiber, projected to the sky. [arcsec]

    efficiency : numpy array
        Spectrograph+telescope throughput for each wavelength bin
        from fig 16 in Kelz+2005. WARNING: Assuming this also includes 
        the Quantum Efficiency. 

    lpmm :  astropy.unit.Quantity
        grating groove frequency [grooves/mm]
        This is a variable for now, until instrument is finalized

    ldisp : astroy.unit.Quantity
        dispersion in Angstrom per pixel [AA/pixel]. 
        Currently set to a constant at all wavelenghts

    dark_current :  astropy.unit.Quantity
        Dark current [electron/s]

    ron : astropy.unit.Quantity
        Readout noise [electron/s]
    
    gain : astroy.unit.Quantity
        Gain of the CCD. Assumed=1 and therefore not explicitly used.

    pixperfib : float
        Number of pixels per fiber. 
        Currently hardcoded here
    
    minwave : astropy.unit.Quantity
        Minimum wavelength of GAMAICA [Angstrom]
    
     maxwave: astropy.unit.Quantity
        maximum wavelength of GAMAICA [Angstrom]
    
    wavelength_perpix: astropy.unit.Quantity array
        wavelength per pixel array of GAMAICA, built with the provided ldisp value [Angstrom]


    wavelength_perfib: astropy.unit.Quantity array
        wavelength per fiber array of GAMAICA, built with the ldisp and pixperfib values  [Angstrom]
        Not actually used yet

    efficiency_perpix: astropy.unit.Quantity array
        Total efficiency per pixel of GAMAICA

    efficiency_perfib: astropy.unit.Quantity array
        Total efficiency per fiber of GAMAICA
        Not actually used yet
    '''
    
    def __init__(self,lpmm):
          
        
        #self.dwave = defaults.lpmm2ldisp[lpmm] * u.Unit('fiber')   # delta wavelength per fiber
        self.ldisp = defaults.lpmm2ldisp[lpmm]   # dispersion in Angstroms per pixel
        self.pixperfib = 4 # pixels per fiber
        self.dark_current = defaults.dark            # CCD dark current per pixel [electrons /s]
        self.ron = defaults.ron                      # CCD read-out noise [electrons]
        self.A_tel = defaults.A_tel                  # effective light collecting area (CAHA3.5m:  8.45 m**2)
        self.fib_diam = defaults.fib_diam            # projected fiber diameter [arcsec]
        self.fib_area = defaults.fib_area            # projected fiber area [arcsec2]
        self.gain = defaults.gain                    # CCD gain
        
        tab = QTable.read(defaults.datadir / 'GAMAICA_ETC_throughput.fits')
        
        ## rebin to match dispersion
        self.minwave = defaults.minwave
        self.maxwave = defaults.maxwave

        self.wavelength_perpix = np.arange(self.minwave.value,self.maxwave.value+self.ldisp.value,self.ldisp.value)*self.minwave.unit 
        self.wavelength_perfib = np.arange(self.minwave.value,self.maxwave.value+self.ldisp.value*self.pixperfib,self.ldisp.value*self.pixperfib) *self.minwave.unit 
        
        
        self.efficiency_perpix = rebin_1d_trans_box(tab['transmission'], tab['wavelength'].to('Angstrom').value, self.wavelength_perpix.value) * u.electron/u.photon
        self.efficiency_perfib = rebin_1d_trans_box(tab['transmission'], tab['wavelength'].to('Angstrom').value, self.wavelength_perfib.value) * u.electron/u.photon

        
        
#fiber_unit = u.def_unit(['fiber','4 *u.pixel'],doc='one fiber is 4 pixels')
#u.add_enabled_units(fiber_unit)
