import numpy as np
from astropy.coordinates import SkyCoord, ICRS, BaseCoordinateFrame

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