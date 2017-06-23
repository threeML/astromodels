import numpy as np

def vincenty(lon0, lat0, a1, s):
  """
  Returns the coordinates of a new point that is a given angular distance s away from a starting point (lon0, lat0) at bearing (angle from north) a1), to within a given precision

  Note that this calculation is a simplified version of the full vincenty problem, which solves for the coordinates on the surface on an arbitrary ellipsoid. Here we only care about the surface of a sphere.

  Note: All parameters are assumed to be given in DEGREES
  :param lon0: float, longitude of starting point
  :param lat0: float, latitude of starting point
  :param a1: float, bearing to second point, i.e. angle between due north and line connecting 2 points
  :param s: float, angular distance between the two points
  :return: coordinates of second point in degrees
  """

  lon0 = np.deg2rad(lon0)
  lat0 = np.deg2rad(lat0)
  a1   = np.deg2rad(a1)
  s    = np.deg2rad(s)

  sina = np.cos(lat0) * np.sin(a1)

  num1 = np.sin(lat0)*np.cos(s) + np.cos(lat0)*np.sin(s)*np.cos(a1)
  den1 = np.sqrt(sina**2 + (np.sin(lat0)*np.sin(s) - np.cos(lat0)*np.cos(a1))**2)
  lat  = np.rad2deg(np.arctan2(num1, den1))

  num2 = np.sin(s)*np.sin(a1)
  den2 = np.cos(lat0)*np.cos(s) - np.sin(lat0)*np.sin(s)*np.cos(a1)
  L    = np.arctan2(num2, den2)
  lon  = np.rad2deg(lon0 + L)

  return lon, lat
