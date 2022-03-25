from __future__ import division

import hashlib

import astropy.units as u
import numpy as np
from astropy import wcs
from astropy.coordinates import ICRS, BaseCoordinateFrame, SkyCoord
from astropy.io import fits
from past.utils import old_div

from astromodels.functions.function import Function2D, FunctionMeta
from astromodels.utils.angular_distance import angular_distance
from astromodels.utils.logging import setup_logger
from astromodels.utils.vincenty import vincenty

log = setup_logger(__name__)

class Latitude_galactic_diffuse(Function2D, metaclass=FunctionMeta):
    r"""
        description :

            A Gaussian distribution in Galactic latitude around the Galactic plane

        latex : $ K \exp{\left( \frac{-b^2}{2 \sigma_b^2} \right)} $

        parameters :

            K :

                desc : normalization
                initial value : 1

            sigma_b :

                desc : Sigma for
                initial value : 1

            l_min :

                desc : min Longtitude
                initial value : 10

            l_max :

                desc : max Longtitude
                initial value : 30

        """

    # This is optional, and it is only needed if we need more setup after the
    # constructor provided by the meta class

    def _setup(self):

        self._frame = ICRS()

    def set_frame(self, new_frame):
        """
        Set a new frame for the coordinates (the default is ICRS J2000)

        :param new_frame: a coordinate frame from astropy
        :return: (none)
        """
        assert isinstance(new_frame, BaseCoordinateFrame)

        self._frame = new_frame

    def _set_units(self, x_unit, y_unit, z_unit):

        self.K.unit = z_unit
        self.sigma_b.unit = x_unit
        self.l_min.unit = y_unit
        self.l_max.unit = y_unit

    def evaluate(self, x, y, K, sigma_b, l_min, l_max):

        # We assume x and y are R.A. and Dec
        _coord = SkyCoord(ra=x, dec=y, frame=self._frame, unit="deg")

        b = _coord.transform_to('galactic').b.value
        l = _coord.transform_to('galactic').l.value

        return K * np.exp(old_div(-b ** 2, (2 * sigma_b ** 2))) * np.logical_or(np.logical_and(l > l_min, l < l_max),np.logical_and(l_min > l_max, np.logical_or(l > l_min, l < l_max)))

    def get_boundaries(self):

        max_b = self.sigma_b.max_value
        l_min = self.l_min.value
        l_max = self.l_max.value

        _coord = SkyCoord(l=[l_min, l_min, l_max, l_max], b=[max_b * -2., max_b * 2., max_b * 2., max_b * -2.], frame="galactic", unit="deg")

        # no dealing with 0 360 overflow
        min_lat = min(_coord.transform_to("icrs").dec.value)
        max_lat = max(_coord.transform_to("icrs").dec.value)
        min_lon = min(_coord.transform_to("icrs").ra.value)
        max_lon = max(_coord.transform_to("icrs").ra.value)

        return (min_lon, max_lon), (min_lat, max_lat)
        
    def get_total_spatial_integral(self, z=None):  
        """
        Returns the total integral (for 2D functions) or the integral over the spatial components (for 3D functions).
        needs to be implemented in subclasses.

        :return: an array of values of the integral (same dimension as z).
        """

        dL= self.l_max.value-self.l_min.value if self.l_max.value > self.l_min.value else 360 + self.l_max.value - self.l_max.value

        #integral -inf to inf exp(-b**2 / 2*sigma_b**2 ) db = sqrt(2pi)*sigma_b 
        #Note that K refers to the peak diffuse flux (at b = 0) per square degree.
        integral = np.sqrt( 2*np.pi ) * self.sigma_b.value * self.K.value * dL 

        if isinstance( z, u.Quantity):
            z = z.value
        return integral * np.power( 180. / np.pi, -2 ) * np.ones_like( z )
        

class Gaussian_on_sphere(Function2D, metaclass=FunctionMeta):
    r"""
        description :

            A bidimensional Gaussian function on a sphere (in spherical coordinates)

        latex : $$ f(\vec{x}) = \left(\frac{180^\circ}{\pi}\right)^2 \frac{1}{2\pi \sqrt{\det{\Sigma}}} \, {\rm exp}\left( -\frac{1}{2} (\vec{x}-\vec{x}_0)^\intercal \cdot \Sigma^{-1}\cdot (\vec{x}-\vec{x}_0)\right) \\ \vec{x}_0 = ({\rm RA}_0,{\rm Dec}_0)\\ \Lambda = \left( \begin{array}{cc} \sigma^2 & 0 \\ 0 & \sigma^2 (1-e^2) \end{array}\right) \\ U = \left( \begin{array}{cc} \cos \theta & -\sin \theta \\ \sin \theta & cos \theta \end{array}\right) \\\Sigma = U\Lambda U^\intercal $$

        parameters :

            lon0 :

                desc : Longitude of the center of the source
                initial value : 0.0
                min : 0.0
                max : 360.0

            lat0 :

                desc : Latitude of the center of the source
                initial value : 0.0
                min : -90.0
                max : 90.0

            sigma :

                desc : Standard deviation of the Gaussian distribution
                initial value : 10
                min : 0
                max : 20

        """

    def _set_units(self, x_unit, y_unit, z_unit):

        # lon0 and lat0 and rdiff have most probably all units of degrees. However,
        # let's set them up here just to save for the possibility of using the
        # formula with other units (although it is probably never going to happen)

        self.lon0.unit = x_unit
        self.lat0.unit = y_unit
        self.sigma.unit = x_unit

    def evaluate(self, x, y, lon0, lat0, sigma):

        lon, lat = x,y

        angsep = angular_distance(lon0, lat0, lon, lat)

        s2 = sigma**2

        return (old_div(180, np.pi))**2 * 1 / (2.0 * np.pi * s2) * np.exp(-0.5 * angsep**2/s2)

    def get_boundaries(self):

        # Truncate the gaussian at 2 times the max of sigma allowed

        max_sigma = self.sigma.max_value

        min_lat = max(-90., self.lat0.value - 2 * max_sigma)
        max_lat = min(90., self.lat0.value + 2 * max_sigma)

        max_abs_lat = max(np.absolute(min_lat), np.absolute(max_lat))

        if max_abs_lat > 89. or 2 * max_sigma / np.cos(max_abs_lat * np.pi / 180.) >= 180.:

            min_lon = 0.
            max_lon = 360.

        else:

            min_lon = self.lon0.value - 2 * max_sigma / np.cos(max_abs_lat * np.pi / 180.)
            max_lon = self.lon0.value + 2 * max_sigma / np.cos(max_abs_lat * np.pi / 180.)

            if min_lon < 0.:

                min_lon = min_lon + 360.

            elif max_lon > 360.:

                max_lon = max_lon - 360.

        return (min_lon, max_lon), (min_lat, max_lat)

    def get_total_spatial_integral(self, z=None):  
        """
        Returns the total integral (for 2D functions) or the integral over the spatial components (for 3D functions).
        needs to be implemented in subclasses.

        :return: an array of values of the integral (same dimension as z).
        """
        
        if isinstance( z, u.Quantity):
            z = z.value
        return np.ones_like( z )
        


class Asymm_Gaussian_on_sphere(Function2D, metaclass=FunctionMeta):
    r"""
        description :

            A bidimensional Gaussian function on a sphere (in spherical coordinates)

            see https://en.wikipedia.org/wiki/Gaussian_function#Two-dimensional_Gaussian_function
            
        parameters :

            lon0 :

                desc : Longitude of the center of the source
                initial value : 0.0
                min : 0.0
                max : 360.0

            lat0 :

                desc : Latitude of the center of the source
                initial value : 0.0
                min : -90.0
                max : 90.0

            a :

                desc : Standard deviation of the Gaussian distribution (major axis)
                initial value : 10
                min : 0
                max : 20

            e :

                desc : Excentricity of Gaussian ellipse
                initial value : 0.9
                min : 0
                max : 1
                
            theta :

                desc : inclination of major axis to a line of constant latitude
                initial value : 10.
                min : -90.0
                max : 90.0

        """

    def _set_units(self, x_unit, y_unit, z_unit):

        # lon0 and lat0 and a have most probably all units of degrees. However,
        # let's set them up here just to save for the possibility of using the
        # formula with other units (although it is probably never going to happen)

        self.lon0.unit = x_unit
        self.lat0.unit = y_unit
        self.a.unit = x_unit
        self.e.unit = u.dimensionless_unscaled
        self.theta.unit = u.degree

    def evaluate(self, x, y, lon0, lat0, a, e, theta):

        lon, lat = x,y

        b = a * np.sqrt(1. - e**2)
        
        dX = np.atleast_1d( angular_distance( lon0, lat0, lon, lat0) )
        dY = np.atleast_1d( angular_distance( lon0, lat0, lon0, lat) )

        dlon = lon - lon0
        if isinstance( dlon, u.Quantity):
            dlon = (dlon.to(u.degree)).value
        
        idx=np.logical_and( np.logical_or( dlon < 0, dlon > 180), np.logical_or( dlon>-180, dlon < -360) )
        dX[idx] = -dX[idx]
        
        idx = lat < lat0 
        dY[idx]=-dY[idx]

        if isinstance( theta, u.Quantity ):
            phi = (theta.to(u.degree)).value + 90.0
        else:
            phi = theta + 90.

        cos2_phi = np.power( np.cos( phi * np.pi/180.), 2)
        sin2_phi = np.power( np.sin( phi * np.pi/180.), 2)
        
        sin_2phi = np.sin( 2. * phi * np.pi/180.)
        
        A = old_div(cos2_phi, (2.*b**2)) + old_div(sin2_phi, (2.*a**2))

        B = old_div(- sin_2phi, (4.*b**2)) + old_div(sin_2phi, (4.*a**2))

        C = old_div(sin2_phi, (2.*b**2)) + old_div(cos2_phi, (2.*a**2))

        E = -A*np.power(dX, 2) + 2.*B*dX*dY - C*np.power(dY, 2)

        return np.power(old_div(180, np.pi), 2) * 1. / (2 * np.pi * a * b) * np.exp( E )
        
    def get_boundaries(self):

        # Truncate the gaussian at 2 times the max of sigma allowed

        min_lat = max(-90., self.lat0.value - 2 * self.a.max_value)
        max_lat = min(90., self.lat0.value + 2 * self.a.max_value)

        max_abs_lat = max(np.absolute(min_lat), np.absolute(max_lat))

        if max_abs_lat > 89. or 2 * self.a.max_value / np.cos(max_abs_lat * np.pi / 180.) >= 180.:

            min_lon = 0.
            max_lon = 360.

        else:

            min_lon = self.lon0.value - 2 * self.a.max_value / np.cos(max_abs_lat * np.pi / 180.)
            max_lon = self.lon0.value + 2 * self.a.max_value / np.cos(max_abs_lat * np.pi / 180.)

            if min_lon < 0.:

                min_lon = min_lon + 360.

            elif max_lon > 360.:

                max_lon = max_lon - 360.

        return (min_lon, max_lon), (min_lat, max_lat)

    def get_total_spatial_integral(self, z=None):  
        """
        Returns the total integral (for 2D functions) or the integral over the spatial components (for 3D functions).
        needs to be implemented in subclasses.

        :return: an array of values of the integral (same dimension as z).
        """

        if isinstance( z, u.Quantity):
            z = z.value
        return np.ones_like( z )


class Disk_on_sphere(Function2D, metaclass=FunctionMeta):
    r"""
        description :

            A bidimensional disk/tophat function on a sphere (in spherical coordinates)

        latex : $$ f(\vec{x}) = \left(\frac{180}{\pi}\right)^2 \frac{1}{\pi~({\rm radius})^2} ~\left\{\begin{matrix} 1 & {\rm if}& {\rm | \vec{x} - \vec{x}_0| \le {\rm radius}} \\ 0 & {\rm if}& {\rm | \vec{x} - \vec{x}_0| > {\rm radius}} \end{matrix}\right. $$

        parameters :

            lon0 :

                desc : Longitude of the center of the source
                initial value : 0.0
                min : 0.0
                max : 360.0

            lat0 :

                desc : Latitude of the center of the source
                initial value : 0.0
                min : -90.0
                max : 90.0

            radius :

                desc : Radius of the disk
                initial value : 15
                min : 0
                max : 20

        """

    def _set_units(self, x_unit, y_unit, z_unit):

        # lon0 and lat0 and rdiff have most probably all units of degrees. However,
        # let's set them up here just to save for the possibility of using the
        # formula with other units (although it is probably never going to happen)

        self.lon0.unit = x_unit
        self.lat0.unit = y_unit
        self.radius.unit = x_unit

    def evaluate(self, x, y, lon0, lat0, radius):

        lon, lat = x,y

        angsep = angular_distance(lon0, lat0, lon, lat)

        return np.power(old_div(180, np.pi), 2) * 1. / (np.pi * radius ** 2) * (angsep <= radius)

    def get_boundaries(self):

        # Truncate the disk at 2 times the max of radius allowed

        max_radius = self.radius.max_value

        min_lat = max(-90., self.lat0.value - 2 * max_radius)
        max_lat = min(90., self.lat0.value + 2 * max_radius)

        max_abs_lat = max(np.absolute(min_lat), np.absolute(max_lat))

        if max_abs_lat > 89. or 2 * max_radius / np.cos(max_abs_lat * np.pi / 180.) >= 180.:

            min_lon = 0.
            max_lon = 360.

        else:

            min_lon = self.lon0.value - 2 * max_radius / np.cos(max_abs_lat * np.pi / 180.)
            max_lon = self.lon0.value + 2 * max_radius / np.cos(max_abs_lat * np.pi / 180.)

            if min_lon < 0.:

                min_lon = min_lon + 360.

            elif max_lon > 360.:

                max_lon = max_lon - 360.

        return (min_lon, max_lon), (min_lat, max_lat)

    def get_total_spatial_integral(self, z=None):  
        """
        Returns the total integral (for 2D functions) or the integral over the spatial components (for 3D functions).
        needs to be implemented in subclasses.

        :return: an array of values of the integral (same dimension as z).
        """

        if isinstance( z, u.Quantity):
            z = z.value
        return np.ones_like( z )


class Ellipse_on_sphere(Function2D, metaclass=FunctionMeta):
    r"""
        description :

            An ellipse function on a sphere (in spherical coordinates)

        latex : $$ f(\vec{x}) = \left(\frac{180}{\pi}\right)^2 \frac{1}{\pi~ a b} ~\left\{\begin{matrix} 1 & {\rm if}& {\rm | \vec{x} - \vec{x}_{f1}| + | \vec{x} - \vec{x}_{f2}| \le {\rm 2a}} \\ 0 & {\rm if}& {\rm | \vec{x} - \vec{x}_{f1}| + | \vec{x} - \vec{x}_{f2}| > {\rm 2a}} \end{matrix}\right. $$

        parameters :

            lon0 :

                desc : Longitude of the center of the source
                initial value : 0.0
                min : 0.0
                max : 360.0

            lat0 :

                desc : Latitude of the center of the source
                initial value : 0.0
                min : -90.0
                max : 90.0

            a :

                desc : semimajor axis of the ellipse
                initial value : 15.
                min : 0
                max : 20
                
            e :

                desc : eccentricity of ellipse 
                initial value : 0.9
                min : 0
                max : 1

            theta :

                desc : inclination of semimajoraxis to a line of constant latitude
                initial value : 0.0
                min : -90.0
                max : 90.0
        """
    
    lon1 = None
    lat1 = None
    lon2 = None
    lat2 = None
    focal_pts = False

    def _set_units(self, x_unit, y_unit, z_unit):

        # lon0 and lat0 have most probably all units of degrees.
        # However, let's set them up here just to save for the possibility of
        # using the formula with other units (although it is probably never
        # going to happen)

        self.lon0.unit = x_unit
        self.lat0.unit = y_unit
        self.a.unit = x_unit
        # eccentricity is dimensionless
        self.e.unit = u.dimensionless_unscaled 
        self.theta.unit = u.degree
        
    def calc_focal_pts(self, lon0, lat0, a, b, theta):
        # focal distance
        f = np.sqrt(a**2 - b**2)

        if isinstance( theta, u.Quantity):
          bearing = 90. - (theta.to(u.degree)).value
        else:
          bearing = 90. - theta
          
        # coordinates of focal points
        lon1, lat1 = vincenty(lon0, lat0, bearing, f)
        lon2, lat2 = vincenty(lon0, lat0, bearing + 180., f)

        return lon1, lat1, lon2, lat2

    def evaluate(self, x, y, lon0, lat0, a, e, theta):
        
        b = a * np.sqrt(1. - e**2)
        
        # calculate focal points
        
        self.lon1, self.lat1, self.lon2, self.lat2 = self.calc_focal_pts(lon0, lat0, a, b, theta)
        self.focal_pts = True
        
        # lon/lat of point in question
        lon, lat = x, y
        
        # sum of geodesic distances to focii (should be <= 2a to be in ellipse)
        angsep1 = angular_distance(self.lon1, self.lat1, lon, lat)
        angsep2 = angular_distance(self.lon2, self.lat2, lon, lat)
        angsep  = angsep1 + angsep2
        
        return np.power(old_div(180, np.pi), 2) * 1. / (np.pi * a * b) * (angsep <= 2*a)

    def get_boundaries(self):

        # Truncate the ellipse at 2 times the max of semimajor axis allowed

        max_radius = self.a.max_value

        min_lat = max(-90., self.lat0.value - 2 * max_radius)
        max_lat = min(90., self.lat0.value + 2 * max_radius)

        max_abs_lat = max(np.absolute(min_lat), np.absolute(max_lat))

        if max_abs_lat > 89. or 2 * max_radius / np.cos(max_abs_lat * np.pi / 180.) >= 180.:

            min_lon = 0.
            max_lon = 360.

        else:

            min_lon = self.lon0.value - 2 * max_radius / np.cos(max_abs_lat * np.pi / 180.)
            max_lon = self.lon0.value + 2 * max_radius / np.cos(max_abs_lat * np.pi / 180.)

            if min_lon < 0.:

                min_lon = min_lon + 360.

            elif max_lon > 360.:

                max_lon = max_lon - 360.

        return (min_lon, max_lon), (min_lat, max_lat)

    def get_total_spatial_integral(self, z=None):  
        """
        Returns the total integral (for 2D functions) or the integral over the spatial components (for 3D functions).
        needs to be implemented in subclasses.

        :return: an array of values of the integral (same dimension as z).
        """

        if isinstance( z, u.Quantity):
            z = z.value
        return np.ones_like( z )


class SpatialTemplate_2D(Function2D, metaclass=FunctionMeta):
    r"""
        description :
        
            User input Spatial Template.  Expected to be normalized to 1/sr
        
        latex : $ hi $
        
        parameters :
        
            K :
        
                desc : normalization
                initial value : 1
                fix : yes
            hash :
                
                desc: hash of model map [needed for memoization]
                initial value: 1
                fix: yes
            ihdu:
                desc: header unit index of fits file
                initial value: 0
                fix: True 
                min: 0
        
        properties:
            fits_file:
                desc: fits file to load
                defer: True
                function: _load_file
            frame:
                desc: coordinate frame
                initial value: icrs
                allowed values:
                    - icrs
                    - galactic
                    - fk5
                    - fk4
                    - fk4_no_e
    """
    
    def _set_units(self, x_unit, y_unit, z_unit):
        
        self.K.unit = z_unit
    
    # This is optional, and it is only needed if we need more setup after the
    # constructor provided by the meta class
    
    # def _setup(self):
        
        # self._frame = "icrs"
        # self._fitsfile = None
        # self._map = None
    
    def _load_file(self):

        self._fitsfile=self.fits_file.value
        
        with fits.open(self._fitsfile) as f:
    
            self._wcs = wcs.WCS( header = f[int(self.ihdu.value)].header )
            self._map = f[int(self.ihdu.value)].data
              
            self._nX = f[int(self.ihdu.value)].header['NAXIS1']
            self._nY = f[int(self.ihdu.value)].header['NAXIS2']

            #note: map coordinates are switched compared to header. NAXIS1 is coordinate 1, not 0. 
            #see http://docs.astropy.org/en/stable/io/fits/#working-with-image-data
            assert self._map.shape[1] == self._nX, "NAXIS1 = %d in fits header, but %d in map" % (self._nX, self._map.shape[1])
            assert self._map.shape[0] == self._nY, "NAXIS2 = %d in fits header, but %d in map" % (self._nY, self._map.shape[0])
            
            #test if the map is normalized as expected
            area = wcs.utils.proj_plane_pixel_area( self._wcs )
            dOmega = (area*u.deg*u.deg).to(u.sr).value
            total = self._map.sum() * dOmega

            if not np.isclose( total, 1,  rtol=1e-2):
                log.warning("2D template read from {} is normalized to {} (expected: 1)".format(self._fitsfile, total) )
            
            #hash sum uniquely identifying the template function (defined by its 2D map array and coordinate system)
            #this is needed so that the memoization won't confuse different SpatialTemplate_2D objects.
            h = hashlib.sha224()
            h.update( self._map)
            h.update( repr(self._wcs).encode('utf-8') )
            self.hash = int(h.hexdigest(), 16)
            

    # def to_dict(self, minimal=False):

    #      data = super(Function2D, self).to_dict(minimal)

    #      if not minimal:
         
    #         data['extra_setup'] = {"_fitsfile": self._fitsfile, "_frame": self._frame }
  
    #      return data
        
    
    # def set_frame(self, new_frame):
    #     """
    #         Set a new frame for the coordinates (the default is ICRS J2000)
            
    #         :param new_frame: a coordinate frame from astropy
    #         :return: (none)
    #         """
    #     assert new_frame.lower() in ['icrs', 'galactic', 'fk5', 'fk4', 'fk4_no_e' ]
                
    #     self._frame = new_frame
    
    def evaluate(self, x, y, K, hash, ihdu):
                  
        # We assume x and y are R.A. and Dec
        coord = SkyCoord(ra=x, dec=y, frame=self.frame.value, unit="deg")
        
        #transform input coordinates to pixel coordinates; 
        #SkyCoord takes care of necessary coordinate frame transformations.
        Xpix, Ypix = coord.to_pixel(self._wcs)
        
        Xpix = np.atleast_1d(Xpix.astype(int))
        Ypix = np.atleast_1d(Ypix.astype(int))
        
        # find pixels that are in the template ROI, otherwise return zero
        #iz = np.where((Xpix<self._nX) & (Xpix>=0) & (Ypix<self._nY) & (Ypix>=0))[0]
        iz = (Xpix<self._nX) & (Xpix>=0) & (Ypix<self._nY) & (Ypix>=0)
        out = np.zeros_like(Xpix).astype(float)
        out[iz] = self._map[Ypix[iz], Xpix[iz]]
        
        return np.multiply(K, out)

    def get_boundaries(self):
    
        # if self._map is None:
            
        #     self.load_file(self._fitsfile)
          
        #We use the max/min RA/Dec of the image corners to define the boundaries.
        #Use the 'outside' of the pixel corners, i.e. from pixel 0 to nX in 0-indexed accounting.
    
        Xcorners = np.array( [0, 0,        self._nX, self._nX] )
        Ycorners = np.array( [0, self._nY, 0,        self._nY] )
        
        corners = SkyCoord.from_pixel( Xcorners, Ycorners, wcs=self._wcs, origin = 0).transform_to(self.frame.value)  
     
        min_lon = min(corners.ra.degree)
        max_lon = max(corners.ra.degree)
        
        min_lat = min(corners.dec.degree)
        max_lat = max(corners.dec.degree)
        
        return (min_lon, max_lon), (min_lat, max_lat)



    def get_total_spatial_integral(self, z=None):  
        """
        Returns the total integral (for 2D functions) or the integral over the spatial components (for 3D functions).
        needs to be implemented in subclasses.

        :return: an array of values of the integral (same dimension as z).
        """

        if isinstance( z, u.Quantity):
            z = z.value
        return np.multiply(self.K.value,np.ones_like( z ) )
        
class Power_law_on_sphere(Function2D, metaclass=FunctionMeta):        

    r"""
        description :

            A power law function on a sphere (in spherical coordinates)

        latex : $$ f(\vec{x}) = \left(\frac{180}{\pi}\right)^{-1.*index}  \left\{\begin{matrix} 0.05^{index} & {\rm if} & |\vec{x}-\vec{x}_0| \le 0.05\\ |\vec{x}-\vec{x}_0|^{index} & {\rm if} & 0.05 < |\vec{x}-\vec{x}_0| \le maxr \\ 0 & {\rm if} & |\vec{x}-\vec{x}_0|>maxr\end{matrix}\right. $$

        parameters :

            lon0 :

                desc : Longitude of the center of the source
                initial value : 0.0
                min : 0.0
                max : 360.0

            lat0 :

                desc : Latitude of the center of the source
                initial value : 0.0
                min : -90.0
                max : 90.0

            index :

                desc : power law index
                initial value : -2.0
                min : -5.0
                max : -1.0

            maxr :

                desc : max radius
                initial value : 20.
                fix : yes

            minr :

                desc : radius below which the PL is approximated as a constant
                initial value : 0.05
                fix : yes

        """

    def _set_units(self, x_unit, y_unit, z_unit):

        # lon0 and lat0 and rdiff have most probably all units of degrees. However,
        # let's set them up here just to save for the possibility of using the
        # formula with other units (although it is probably never going to happen)

        self.lon0.unit = x_unit
        self.lat0.unit = y_unit
        self.index.unit = u.dimensionless_unscaled
        self.maxr.unit = x_unit
        self.minr.unit = x_unit

    def evaluate(self, x, y, lon0, lat0, index, maxr, minr):

        lon, lat = x,y

        angsep = angular_distance(lon0, lat0, lon, lat)

        if maxr <= minr:
            norm = np.power(np.pi / 180., 2.+index) * np.pi * maxr**2 * minr**index
        elif self.index.value == -2.:
            norm = np.pi * (1.0 + 2. * np.log(maxr / minr) )
        else:
            norm = np.power(minr * np.pi / 180., 2.+index) * np.pi + 2. * np.pi / (2.+index) * (np.power(maxr * np.pi / 180., index+2.) - np.power(minr * np.pi / 180., index+2.))

        value = np.less_equal(angsep,maxr) * np.power(np.pi / 180., index) * np.power(np.add(np.multiply(angsep, np.greater(angsep, minr)), np.multiply(minr, np.less_equal(angsep, minr))), index)

        return value / norm

    def get_boundaries(self):

        return ((self.lon0.value - self.maxr.value), (self.lon0.value + self.maxr.value)), ((self.lat0.value - self.maxr.value), (self.lat0.value + self.maxr.value))

    def get_total_spatial_integral(self, z=None):  
        """
        Returns the total integral (for 2D functions) or the integral over the spatial components (for 3D functions).
        needs to be implemented in subclasses.

        :return: an array of values of the integral (same dimension as z).
        """
        
        if isinstance( z, u.Quantity):
            z = z.value
        return np.ones_like( z )
 

# class FunctionIntegrator(Function2D):
#     r"""
#         description :
#
#             Returns the average of the integrand function (a 1-d function) over the interval x-y. The integrand is set
#             using the .integrand property, like in:
#
#             > G = FunctionIntegrator()
#             > G.integrand = Powerlaw()
#
#         latex : $$ G(x,y) = \frac{\int_{x}^{y}~f(x)~dx}{y-x}$$
#
#         parameters :
#
#             s :
#
#                 desc : if s=0, then the integral will *not* be normalized by (y-x), otherwise (default) it will.
#                 initial value : 1
#                 fix : yes
#         """
#
#     __metaclass__ = FunctionMeta
#
#     def _set_units(self, x_unit, y_unit, z_unit):
#
#         # lon0 and lat0 and rdiff have most probably all units of degrees. However,
#         # let's set them up here just to save for the possibility of using the
#         # formula with other units (although it is probably never going to happen)
#
#         self.s = u.dimensionless_unscaled
#
#     def evaluate(self, x, y, s):
#
#         assert y-x >= 0, "Cannot obtain the integral if the integration interval is zero or negative!"
#
#         integral = self._integrand.integral(x, y)
#
#         if s==0:
#
#             return integral
#
#         else:
#
#             return integral / (y-x)
#
#
#     def get_boundaries(self):
#
#         return (-np.inf, +np.inf), (-np.inf, +np.inf)
#
#     def _set_integrand(self, function):
#
#         self._integrand = function
#
#     def _get_integrand(self):
#
#         return self._integrand
#
#     integrand = property(_get_integrand, _set_integrand,
#                          doc="""Get/set the integrand""")
#
#
#     def to_dict(self, minimal=False):
#
#         data = super(Function2D, self).to_dict(minimal)
#
#         if not minimal:
#             data['extra_setup'] = {'integrand': self.integrand.path}
#
#         return data
