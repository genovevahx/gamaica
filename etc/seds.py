from .stddata import seddir
from astropy.table import QTable
import pysynphot as S
import glob
import astropy.units as u

class SEDtemplate:
    '''Spectral energy distribution of the target of observation

    Initialization happens using one of the templates in ./seds/ directory
    ToDo: User-provided SED is coded, but not tested

    The target SED must be explicitly scaled to the desired magnitude, see example

    Parameters
    ----------
    name : str
        Name of the template. Must be exactly like the file names in ./seds/ dir but without the '.fits'
        For now, these are available:

        Filename            - `name`
        'Kinney_starb1.fits'- Kinney_starb1
        'orion.fits'        - orion
        'Pickles_A0V.fits'  - Pickles_A0V
        'Pickles_O5V.fits'  - Pickles_O5V
        'pn.fits'           - pn
        'qso.fits'          - qso

    Examples
    --------

    Create an HII region SED and scale it to Vmag = 20 AB

        >>> from etc.seds import SEDtemplate
        >>> HII_region = SEDtemplate('orion')
        >>> HII_region.scale_to_mag(20,'johnson,v')
       

    '''
    def __init__(self,name):

        self.filter_name = ''
        ## look for template in default directory
        fnames = list(seddir.glob('%s.fits'%name))
        if len(fnames)==1:
            fname = fnames[0]
            tab = QTable.read(fname)
            self.sed = S.ArraySpectrum(tab['wavelength'], tab['flam'],waveunits='angstrom',fluxunits='flam')
        else:
            # if default directory empty, it's a user provided template
            try:
                tab = QTable.read(name,format='ascii')
                self.sed = S.ArraySpectrum(tab['col1']*u.Angstrom, tab['col2']*u.erg/u.s/u.cm**2/u.AA,waveunits='angstrom',fluxunits='flam')
            except:
                raise ValueError(f'File name {name} does not exist in seds directory or in present working directory')
        
        
        

    
    def scale_to_mag(self,mag,filter_):
        '''
        scales the SED template to the given magnitude in the given filter
        mag : scalar
           magnitude of the target, assigned to the SED template [ABmag]
        bandpass : str 
           has to be one of 'johnson,v','johnson,b','cousins,r','cousins,i'

        '''
        import astropy.constants as const

        FILTER = S.ObsBandpass(filter_)
        obs = S.Observation(self.sed, FILTER, binset = self.sed.wave*u.Angstrom, force='taper')
        
        integrated_flam_in_filter =  obs.effstim('flam')
        desired_flam_in_filter  = 10**(-0.4*(mag+48.6)) * const.c.to('AA/s').value / FILTER.pivot()**2
        
        scaling_constant = desired_flam_in_filter / integrated_flam_in_filter
        
        #print('Multiplying SED by ',scaling_constant)
        self.sed = S.ArraySpectrum(wave=self.sed.wave, flux=self.sed.flux * scaling_constant, waveunits=self.sed.waveunits, fluxunits=self.sed.fluxunits)
        self.pivot_wave = FILTER.pivot()
