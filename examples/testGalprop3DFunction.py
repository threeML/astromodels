# code to show how to use the GalPropTemplate_3D() function
# author: Hugo Ayala (hgayala@psu.edu)
# date: Apr 22, 2019

import numpy as np
from threeML import GalPropTemplate_3D

MODELMAP = "testICModel.fits"

# Create  and initialze shape object
shape = GalPropTemplate_3D()

# Load template file. Information that is needed is:
# Filename, RA(Long) low, RA(Long) high,
# Dec(Lat) low, Dec(Lat), high, and bool variable to tell if coordinates
# are in galactic or celestial.
# The galprop map was generated with flux information in units of Mev / (cm^2 s sr)
shape.load_file(MODELMAP, 65, 75, -5, 5, True)

# The templates are cubes that contain spatial, energy and flux information.

# Functions that are called during the fitting
print("Function get boundaries()")
(ramin, ramax), (decmin, decmax) = shape.get_boundaries()
print("RA: {:0.2f},{:0.2f} Dec: {:0.2f},{:0.2f}".format(ramin, ramax, decmin, decmax))

print("\nFunction evaluate()")

x = np.linspace(ramin, ramax, 10)
y = np.linspace(decmin, decmax, 10)
e = np.array([100000, 1000000, 10000000, 100000000, 1000000000])
res = shape.evaluate(x, y, e, 1, 1)
print(res)
