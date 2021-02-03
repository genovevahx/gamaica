from etc.atmosphere import Atmosphere
from etc.seds import SEDtemplate
from etc.instrument import Gamaica
from etc.exposure import Exposure
from etc.functions import save_and_see

## instrument
instrument = Gamaica(1200)

## atmosphere
atm = Atmosphere.caha()
sky_emission = atm.emission(20.5,'johnson,v')
sky_transmission = atm.transmission(airmass=1.2)

# target to observe
HII_region = SEDtemplate('orion')
HII_region.scale_to_mag(20.0,'johnson,v')


# set up exposure
exp = Exposure(instrument)

# add object to observe
exp.add_object(HII_region)

# add sky emission
exp.add_sky_emission(sky_emission)

# add atmospheric transmission (optional)
exp.add_sky_transmission(sky_transmission)

# expose
exptime = 3600. # seconds
result = exp.expose(exptime)

print(result)
# save result, look at plots
save_and_see(result)



