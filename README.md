# gamaica

Prerequisites
* pysynphot v1.0.0
* astropy 1.1 or greater
* numpy 1.9 or greater
* matplotlib (optional) 

To install pysynphot: 

pip install pysynphot

from terminal. You then need to download the library pysynphot uses from here:
http://ssb.stsci.edu/trds/tarfiles/synphot2.tar.gz
and unpack it in a directory, e.g. I put it in

/z/genoveva/extras/pysynphot/cdbs

then you need to add that directory to the path

export PYSYN_CDBS=/z/genoveva/extras/pysynphot/cdbs/

which you can do in your bashrc. 


Pysynphot will be phased-out at some point and these steps will not be necessary. 

Things to note:


    *Sky emission*: dark sky emission at CAHA is observational, and obtained from 
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
    
    
    * A constant fiber fraction (amount of star entering the fiber) is assumed to be 80%. 
    * dispersion is assumed constant with wavelength
    * binning hasn't been implemented
    * Gray and Bright sky conditions are for the time being approximated by simply using the dark sky spectrum at CAHA and scaling it to a user-specified magnitude in a user-selected filter
    * the standard Johnson B,V and Cousins R,I are available
    * more SED templates can be added but currently there are only 
      - Kinney_starb1 (starburst galaxy, E(B-V)=0.1)
      - pn (planetary nebula)
      - Orion (HII region)
      - Pickles A0V and O5V stars
      - qso (quasar)
    
