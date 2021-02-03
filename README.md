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

