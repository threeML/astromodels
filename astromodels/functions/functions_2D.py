import numpy as np
from scipy.integrate import quad
from astropy.coordinates import SkyCoord, ICRS, BaseCoordinateFrame
from astropy.io import fits
import astropy.units as u

from astromodels.functions.function import Function2D, FunctionMeta
from astromodels.utils.angular_distance import angular_distance


class Latitude_galactic_diffuse(Function2D):
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

        """

    __metaclass__ = FunctionMeta

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

    def evaluate(self, x, y, K, sigma_b):

        # We assume x and y are R.A. and Dec
        _coord = SkyCoord(ra=x, dec=y, frame=self._frame, unit="deg")

        b = _coord.transform_to('galactic').b.value

        return K * np.exp(-b ** 2 / (2 * sigma_b ** 2))


class Gaussian_on_sphere(Function2D):
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
                initial value : 0.5
                min : 0
                max : 20

        """

    __metaclass__ = FunctionMeta

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

        return np.power(180 / np.pi, 2) * 1. / (2 * np.pi * sigma ** 2) * np.exp(
            -0.5 * np.power(angsep, 2) / sigma ** 2)

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


class Disk_on_sphere(Function2D):
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
                initial value : 0.5
                min : 0
                max : 20

        """

    __metaclass__ = FunctionMeta

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

        return np.power(180 / np.pi, 2) * 1. / (np.pi * radius ** 2) * (angsep <= radius)

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


class Ellipse_on_sphere(Function2D):
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

                desc : semimajor axis of the disk
                initial value : 0.5
                min : 0
                max : 20
                
            b :

                desc : semiminor axis of the disk
                initial value : 0.5
                min : 0
                max : 20
                
            theta :

                desc : inclination of semimajoraxis to a line of constant latitude
                initial value : 0.0
                min : -90.0
                max : 90.0
        """

    __metaclass__ = FunctionMeta

    def _set_units(self, x_unit, y_unit, z_unit):

        # lon0 and lat0 and rdiff have most probably all units of degrees.
        # However, let's set them up here just to save for the possibility of
        # using the formula with other units (although it is probably never
        # going to happen)

        self.lon0.unit = x_unit
        self.lat0.unit = y_unit
        self.a.unit = x_unit
        self.b.unit = x_unit
        self.theta.unit = x_unit

    def evaluate(self, x, y, lon0, lat0, a, b, theta):

        # lon/lat of point in question
        lon, lat = x,y
        
        # focal distance
        f = np.sqrt(a**2 - b**2)
        
        # focus 1 coordinate and distance from focus 1 to point
        lon1 = lon0 - f*np.cos(theta)
        lat1 = lat0 - f*np.sin(theta)
        angsep1 = angular_distance(lon1, lat1, lon, lat)

        # focus 2 coordinate and distance from focus 2 to point
        lon2 = lon0 + f*np.cos(theta)
        lat2 = lat0 + f*np.sin(theta)
        angsep2 = angular_distance(lon2, lat2, lon, lat)

        # sum of distances to focii (should be <= 2a to be in ellipse)
        angsep = angsep1 + angsep2
        
        return np.power(180 / np.pi, 2) * 1. / (np.pi * a * b) * (angsep <= 2*a)

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


class SpatialTemplate_2D(Function2D):
    r"""
        description :
        
            User input Spatial Template.  Expected to be normalized to 1/sr
        
        latex : $ hi $
        
        parameters :
        
            K :
        
                desc : normalization
                initial value : 1
                fix : yes
        
        """
    
    __metaclass__ = FunctionMeta
    
    def _set_units(self, x_unit, y_unit, z_unit):
        
        self.K.unit = z_unit
    
    # This is optional, and it is only needed if we need more setup after the
    # constructor provided by the meta class
    
    def _setup(self):
        
        self._frame = ICRS()
    
    def load_file(self,fitsfile,ihdu=0):
        
        with fits.open(fitsfile) as f:
            self._refXpix = f[ihdu].header['CRPIX1']
            self._refYpix = f[ihdu].header['CRPIX2']
            self._delXpix = f[ihdu].header['CDELT1']
            self._delYpix = f[ihdu].header['CDELT2']
            self._refX = f[ihdu].header['CRVAL1'] # assumed to be RA
            self._refY = f[ihdu].header['CRVAL2'] # assumed to be DEC
            self._map = f[ihdu].data
            self._nX = f[ihdu].header['NAXIS1']
            self._nY = f[ihdu].header['NAXIS2']
    
    def set_frame(self, new_frame):
        """
            Set a new frame for the coordinates (the default is ICRS J2000)
            
            :param new_frame: a coordinate frame from astropy
            :return: (none)
            """
        assert isinstance(new_frame, BaseCoordinateFrame)
        
        self._frame = new_frame
    
    def evaluate(self, x, y, K):
        
        # We assume x and y are R.A. and Dec
        _coord = SkyCoord(ra=x, dec=y, frame=self._frame, unit="deg")
        
        Xpix = np.add(np.divide(np.subtract(x,self._refX),self._delXpix),self._refXpix)
        Ypix = np.add(np.divide(np.subtract(y,self._refY),self._delYpix),self._refYpix)
        
        Xpix = Xpix.astype(int)
        Ypix = Ypix.astype(int)
        
        # find pixels that are in the template ROI, otherwise return zero
        iz = np.where((Xpix<self._nX) & (Xpix>=0) & (Ypix<self._nY) & (Ypix>=0))[0]
        out = np.zeros((len(x)))
        out[iz] = self._map[Xpix[iz].astype(int),Ypix[iz]]
        
        #pdb.set_trace()
        
        return np.multiply(K,out)

    def get_boundaries(self):
    
        min_ra = (0-np.int(self._refXpix))*self._delXpix + self._refX
        max_ra = ((self._nX-1)-np.int(self._refXpix))*self._delXpix + self._refX
        
        min_dec = (0-np.int(self._refYpix))*self._delYpix + self._refY
        max_dec = ((self._nY-1)-np.int(self._refYpix))*self._delYpix + self._refY
        
        min_lon = min([min_ra,max_ra])
        max_lon = max([min_ra,max_ra])
        
        min_lat = min([min_dec,max_dec])
        max_lat = max([min_dec,max_dec])
        
        
        return (min_lon, max_lon), (min_lat, max_lat)


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
