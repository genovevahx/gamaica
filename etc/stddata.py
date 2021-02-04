import pathlib
import astropy.units as u

datadir = pathlib.Path(__file__).parent / '../data'
icrdir = pathlib.Path(__file__).parent / '../icr'
seddir = pathlib.Path(__file__).parent / '../seds'
dark_sky_Vmag_AB = 20.56490531013599 # dark sky AB mag in Johnson V

lpmm2ldisp = {1200:0.45*u.AA,# dispersion in Angstroms per pixel
              600:0.9*u.AA,
              300:1.8*u.AA}#,


A_tel = 8.49 * u.m**2                       # effective light collecting area (CAHA3.5m:  8.45 m**2)
dark  = (1.0/3600.)*u.electron *u.second**-1# CCD dark current per pixel [electrons/s]
ron   = 5.2 * u.electron	            # CCD read-out noise [electrons]
fib_diam = 1.08 * u.arcsec                  # projected fiber diameter [arcsec]
fib_area = 0.916 * u.arcsec**2              # projected fiber area [arcsec2]
minwave = 4300 * u.AA                       # minimum lambda of GAMAICA
maxwave = 6800 * u.AA                       # maximum lambda of GAMAICA
gain = 1.0 	                            # CCD gain [electrons/ADU]
fiber_fraction = 0.8                        # assume constant fraction of point source entering
                                            # fiber. This is actually wavelength and seeing dependent
