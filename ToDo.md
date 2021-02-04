
**basic ToDos:**
* implement dispersion as a function of wavelength
* maybe separate GAMAICA into explicit Telescope and Instrument classes
* implement fiber fraction loss as a function wavelength and seeing
* implement binning
* add more template SEDs
* phase-out pysynphot by implementing ArraySpectrum class locally
  - add more filter transmission curves
* add test functions
* update __init__.py
* implement gray and bright sky emission correctly (currently just scaling dark sky)

**Comments from Ole Streicher on improving the code:**

I very much like the idea of having and SEDtemplate class. However, the "scale_to_mag()" method there looks dangerous to me, since it manipulates the SED itself. If you ever think of parallelization, this will fail (since one thread may already want to scale it to a new magnitude, while an other thread still uses the old values). I would rather let the function return the scaled spectrum instead.

In this function, you use pysynphot.Filter internally but not externally. If you really want to keep pysynphot (resp. its API), then you could also move that to the used, at least optionally:

def scale_to_mag(self, mag, flt):
    if isinstance(flt, str):
        flt = S.ObsBandPass(flt)  # convenience: give the filter by name
    ...
    return S.ArraySpectrum(...)

I would also not handle ABmag manually, but use astropy.units.ABmag instead.

atmosphere.py: I would consider using SEDtemplate to scale the atmosphere instead of duplicating its code. One of the Python zens is: Don't repeat yourself. Just scaling the dark sky however seems not a correct way to simulate brighter skies; I'd think that scaling f.e. the ESO moon flux would be far more accurate.

instrument.py: Is there a reason to unify telescope and spectrograph? I still think that the separation is useful, as they relate to different steps of the processing.

exposure.py: This seems to change the data flow wrt. 4MOST ETC: I would do

 * create a target spectrum, and an atmosphere
 * pass the target spectrum through atm.transmission
 * pass the target spectrum and atm.emission through the telescope and the spectrograph into exposure.py
 * check exposure.py for different exposure times until one is happy

Having the transmission in exposure.py seems to me the wrong place: once the light is on the exposure, it already went through the atmosphere (... and the telescope, and the spectrograph).

What I did here to simplify everything is to create "Observatory" and "Observation" classes (new file observation.py in my code) which combine everything. In 4MOST it looks like this:

    import astropy.units as u
    from etc import SEDTemplate, QMostObservatory
    observatory = QMostObservatory('hrs')
    obs = observatory(45*u.deg, 1.3*u.arcsec, 'gray')
    pickles = SEDTemplate('Pickles_G0V')
    flux = pickles(15*u.ABmag, 'GAIA2r.G')
    obs.add_object(pickles.wavelength, flux, 'point')
    tbl = obs.expose(300*u.s)

which seems to me a good user interface.

__init__.py: In the meantime, I added there all classes that I want to have visible for the application. I'd recommend you to do the same, as it makes the code simpler.

Generally: I would avoid too much debugging output. Most of this is used for testing -- but testing should be done with proper unit tests. The idea of these is that you don't need to manually check them. As a simple example, in sed.py, you have a printout "Multiplying SED by...", which seems that you are not sure whether the scaling works well. Instead of the printout, test the function with a number of known inputs and magnitudes and check that the scaling is correct. Once you did this, changing the code here f.e. to get rid of the manual "10**(-0.4*(mag+48.6))" is not dangerous anymore, since the tests will make sure that it still works as before.
