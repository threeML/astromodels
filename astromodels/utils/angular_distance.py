from __future__ import division
from past.utils import old_div
import numpy as np


def angular_distance_fast(ra1, dec1, ra2, dec2):
    """
    Compute angular distance using the Haversine formula. Use this one when you know you will never ask for points at
    their antipodes. If this is not the case, use the angular_distance function which is slower, but works also for
    antipodes.

    :param lon1:
    :param lat1:
    :param lon2:
    :param lat2:
    :return:
    """

    lon1 = np.deg2rad(ra1)
    lat1 = np.deg2rad(dec1)
    lon2 = np.deg2rad(ra2)
    lat2 = np.deg2rad(dec2)
    dlon = lon2 - lon1
    dlat = lat2 - lat1

    a = np.sin(old_div(dlat,2.0))**2 + np.cos(lat1) * np.cos(lat2) * np.sin(old_div(dlon,2.0))**2
    c = 2 * np.arcsin(np.sqrt(a))
    return np.rad2deg(c)


def angular_distance(ra1, dec1, ra2, dec2):
    """
    Returns the angular distance between two points, two sets of points, or a set of points and one point.

    :param ra1: array or float, longitude of first point(s)
    :param dec1: array or float, latitude of first point(s)
    :param ra2: array or float, longitude of second point(s)
    :param dec2: array or float, latitude of second point(s)
    :return: angular distance(s) in degrees
    """

    # Vincenty formula, slower than the Haversine formula in some cases, but stable also at antipodes

    lon1 = np.deg2rad(ra1)
    lat1 = np.deg2rad(dec1)
    lon2 = np.deg2rad(ra2)
    lat2 = np.deg2rad(dec2)

    sdlon = np.sin(lon2 - lon1)
    cdlon = np.cos(lon2 - lon1)
    slat1 = np.sin(lat1)
    slat2 = np.sin(lat2)
    clat1 = np.cos(lat1)
    clat2 = np.cos(lat2)

    num1 = clat2 * sdlon
    num2 = clat1 * slat2 - slat1 * clat2 * cdlon
    denominator = slat1 * slat2 + clat1 * clat2 * cdlon

    return np.rad2deg(np.arctan2(np.sqrt(num1 ** 2 + num2 ** 2), denominator))


def spherical_angle( ra0, dec0, ra1, dec1, ra2, dec2 ):
    """
    Returns the spherical angle distance between two sets of great circles defined by (ra0, dec0), (ra1, dec1) and (ra0, dec0), (ra2, dec2)

    :param ra0: array or float, longitude of intersection point(s)
    :param dec0: array or float, latitude of intersection point(s)
    :param ra1: array or float, longitude of first point(s)
    :param dec1: array or float, latitude of first point(s)
    :param ra2: array or float, longitude of second point(s)
    :param dec2: array or float, latitude of second point(s)
    :return: spherical angle in degrees
    """
    
    a = np.deg2rad( angular_distance(ra0, dec0, ra1, dec1))
    b = np.deg2rad( angular_distance(ra0, dec0, ra2, dec2))
    c = np.deg2rad( angular_distance(ra2, dec2, ra1, dec1))
    
    #use the spherical law of cosines: https://en.wikipedia.org/wiki/Spherical_law_of_cosines#Rearrangements
    
    numerator = np.atleast_1d( np.cos(c) - np.cos(a) * np.cos(b) )
    denominator = np.atleast_1d( np.sin(a)*np.sin(b) )
    
    return np.where( denominator == 0 , np.zeros( len(denominator)), np.rad2deg( np.arccos( old_div(numerator,denominator))) )
