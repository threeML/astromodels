# code to show how to create an extended source via uploading a FITs image template
# author: Andrea Albert (aalbert@slac.stanford.edu)
# date: Oct 26, 2016

from threeML import *

# the class SpatialTemplate_2D expects a FITs file that contains a header with the following info: reference pixels (e.g. 'CRPIX1'), pixels step in degrees (e.g. 'CDELT1'), RA and DEC values at reference pixel (e.g. 'CRVAL1')

# initialize shape object
shape = SpatialTemplate_2D()
# load in template file
# by default the extension number is set to zero (ihdu = 0)
shape.load_file("exampleDMtemplate.fits", ihdu=0)

# just for example let's assume a powerlaw spectrum
spectrum = Powerlaw()

source = ExtendedSource("M31", spatial_shape=shape, spectral_shape=spectrum)

# The code assumes the template is normalized to 1 sr.  If it isn't be default then you should set the optional normalization (K) appropriately.  The example template is already normalized to 1 sr so we'll keep K set to 1.  Note K is set to 1 and fixed by default, we include the following commands as an example of how to manipulate K
shape.K = 1.0
shape.K.fix = True

# The following are example commands that get called during fitting

# get the edges of the template
(min_ra, max_ra), (min_dec, max_dec) = shape.get_boundaries()

# return the values at various pixels at locations (x,y).  Note the code assumes x=RA (degrees) and y=DEC(degrees).  Note the code will return a value of 0 is the pixel is outside the template ROI...in this example only the 2nd pixel will have a non-zero value
val = shape.evaluate(x=[1.0, 10.0, 10.0], y=[1.0, 40.0, 89.0], K=1)
