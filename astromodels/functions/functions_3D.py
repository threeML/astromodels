from __future__ import division

import hashlib

import astropy.units as u
import numpy as np
from astropy.coordinates import ICRS, BaseCoordinateFrame, SkyCoord
from astropy.io import fits
from past.utils import old_div
from scipy.interpolate import RegularGridInterpolator

from astromodels.functions.function import Function3D, FunctionMeta
from astromodels.utils.angular_distance import angular_distance_fast


class Continuous_injection_diffusion_ellipse(Function3D, metaclass=FunctionMeta):
    r"""
        description :

            Positron and electrons diffusing away from the accelerator

        latex : $\left(\frac{180^\circ}{\pi}\right)^2 \frac{1.2154}{\sqrt{\pi^3} r_{\rm diff} ({\rm angsep} ({\rm x, y, lon_0, lat_0})+0.06 r_{\rm diff} )} \, {\rm exp}\left(-\frac{{\rm angsep}^2 ({\rm x, y, lon_0, lat_0})}{r_{\rm diff} ^2} \right)$

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

            rdiff0 :

                desc : Projected diffusion radius. The maximum allowed value is used to define the truncation radius.
                initial value : 1.0
                min : 0
                max : 20

            delta :

                desc : index for the diffusion coefficient
                initial value : 0.5
                min : 0.3
                max : 0.6
                fix : yes

            b :

                desc : b field strength in uG
                initial value : 3
                min : 1
                max : 10.
                fix : yes

            piv :

                desc : Pivot for the diffusion radius
                initial value : 2e10
                min : 0
                fix : yes

            piv2 :

                desc : Pivot for converting gamma energy to electron energy (always be 1 TeV)
                initial value : 1e9
                min : 0
                fix : yes

            incl :

                desc : inclination of semimajoraxis to a line of constant latitude
                initial value : 0.0
                min : -90.0
                max : 90.0
                fix : yes

            elongation :

                desc : elongation of the ellipse (b/a)
                initial value : 1.
                min : 0.1
                max : 10.

        """

    def _set_units(self, x_unit, y_unit, z_unit, w_unit):

        # lon0 and lat0 and rdiff have most probably all units of degrees. However,
        # let's set them up here just to save for the possibility of using the
        # formula with other units (although it is probably never going to happen)

        self.lon0.unit = x_unit
        self.lat0.unit = y_unit
        self.rdiff0.unit = x_unit

        # Delta is of course unitless

        self.delta.unit = u.dimensionless_unscaled
        self.b.unit = u.dimensionless_unscaled
        self.incl.unit = x_unit
        self.elongation.unit = u.dimensionless_unscaled

        # Piv has the same unit as energy (which is z)

        self.piv.unit = z_unit
        self.piv2.unit = z_unit

    def evaluate(self, x, y, z, lon0, lat0, rdiff0, delta, b, piv, piv2, incl, elongation):

        lon, lat = x, y
        energy = z

        # energy in kev -> TeV.
        # NOTE: the use of piv2 is necessary to preserve dimensional correctness: the logarithm can only be taken
        # of a dimensionless quantity, so there must be a pivot there.

        e_energy_piv2 = 17. * \
            np.power(old_div(energy, piv2), 0.54 + 0.046 *
                     np.log10(old_div(energy, piv2)))
        e_piv_piv2 = 17. * \
            np.power(old_div(piv, piv2), 0.54 + 0.046 *
                     np.log10(old_div(piv, piv2)))

        try:

            rdiff_a = rdiff0 * np.power(old_div(e_energy_piv2, e_piv_piv2), old_div((delta - 1.), 2.)) * \
                np.sqrt(b * b / 8. / np.pi * 0.624 + 0.26 * np.power(1. + 0.0107 * e_piv_piv2, -1.5)) / \
                np.sqrt(b * b / 8. / np.pi * 0.624 + 0.26 *
                        np.power(1. + 0.0107 * e_energy_piv2, -1.5))

        except ValueError:

            # This happens when using units, because astropy.units fails with the message:
            # "ValueError: Quantities and Units may only be raised to a scalar power"

            # Work around the problem with this loop, which is slow but using units is only for testing purposes or
            # single calls, so it shouldn't matter too much
            rdiff_a = np.array([(rdiff0 * np.power(old_div(e_energy_piv2, e_piv_piv2), x)).value for x in (delta - 1.) / 2. * np.sqrt(b * b / 8. / np.pi * 0.624 + 0.26 * np.power(1. + 0.0107 * e_piv_piv2, -1.5)) /
                                np.sqrt(b * b / 8. / np.pi * 0.624 + 0.26 * np.power(1. + 0.0107 * e_energy_piv2, -1.5))]) * rdiff0.unit

        rdiff_b = rdiff_a * elongation

        pi = np.pi

        angsep = angular_distance_fast(lon, lat, lon0, lat0)
        ang = np.arctan2(lat - lat0, (lon - lon0) *
                         np.cos(lat0 * np.pi / 180.))

        theta = np.arctan2(old_div(np.sin(ang-incl*np.pi/180.),
                                   elongation), np.cos(ang-incl*np.pi/180.))

        rdiffs_a, thetas = np.meshgrid(rdiff_a, theta)
        rdiffs_b, angseps = np.meshgrid(rdiff_b, angsep)

        rdiffs = np.sqrt(rdiffs_a ** 2 * np.cos(thetas) **
                         2 + rdiffs_b ** 2 * np.sin(thetas) ** 2)

        results = np.power(old_div(180.0, pi), 2) * 1.22 / (pi * np.sqrt(pi) * rdiffs_a * np.sqrt(
            elongation) * (angseps + 0.06 * rdiffs)) * np.exp(old_div(-np.power(angseps, 2), rdiffs ** 2))

        return results

    def get_boundaries(self):

        # Truncate the function at the max of rdiff allowed

        maximum_rdiff = self.rdiff0.max_value

        min_latitude = max(-90., self.lat0.value - maximum_rdiff)
        max_latitude = min(90., self.lat0.value + maximum_rdiff)

        max_abs_lat = max(np.absolute(min_latitude), np.absolute(max_latitude))

        if max_abs_lat > 89. or old_div(maximum_rdiff, np.cos(max_abs_lat * np.pi / 180.)) >= 180.:

            min_longitude = 0.
            max_longitude = 360.

        else:

            min_longitude = self.lon0.value - \
                old_div(maximum_rdiff, np.cos(max_abs_lat * np.pi / 180.))
            max_longitude = self.lon0.value + \
                old_div(maximum_rdiff, np.cos(max_abs_lat * np.pi / 180.))

            if min_longitude < 0.:

                min_longitude += 360.

            elif max_longitude > 360.:

                max_longitude -= 360.

        return (min_longitude, max_longitude), (min_latitude, max_latitude)

    def get_total_spatial_integral(self, z=None):
        """
        Returns the total integral (for 2D functions) or the integral over the spatial components (for 3D functions).
        needs to be implemented in subclasses.

        :return: an array of values of the integral (same dimension as z).
        """

        if isinstance(z, u.Quantity):
            z = z.value
        return np.ones_like(z)


class Continuous_injection_diffusion(Function3D, metaclass=FunctionMeta):
    r"""
        description :

            Positron and electrons diffusing away from the accelerator

        latex : $\left(\frac{180^\circ}{\pi}\right)^2 \frac{1.2154}{\sqrt{\pi^3} r_{\rm diff} ({\rm angsep} ({\rm x, y, lon_0, lat_0})+0.06 r_{\rm diff} )} \, {\rm exp}\left(-\frac{{\rm angsep}^2 ({\rm x, y, lon_0, lat_0})}{r_{\rm diff} ^2} \right)$

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

            rdiff0 :

                desc : Projected diffusion radius limited by the cooling time. The maximum allowed value is used to define the truncation radius.
                initial value : 1.0
                min : 0
                max : 20

            rinj :

                desc : Ratio of diffusion radius limited by the injection time over rdiff0. The maximum allowed value is used to define the truncation radius.
                initial value : 100.0
                min : 0
                max : 200
                fix : yes

            delta :

                desc : index for the diffusion coefficient
                initial value : 0.5
                min : 0.3
                max : 0.6
                fix : yes

            b :

                desc : b field strength in uG
                initial value : 3
                min : 1
                max : 10.
                fix : yes

            piv :

                desc : Pivot for the diffusion radius
                initial value : 2e10
                min : 0
                fix : yes

            piv2 :
                desc : Pivot for converting gamma energy to electron energy (always be 1 TeV)
                initial value : 1e9
                min : 0
                fix : yes

        """

    def _set_units(self, x_unit, y_unit, z_unit, w_unit):

        # lon0 and lat0 and rdiff have most probably all units of degrees. However,
        # let's set them up here just to save for the possibility of using the
        # formula with other units (although it is probably never going to happen)

        self.lon0.unit = x_unit
        self.lat0.unit = y_unit
        self.rdiff0.unit = x_unit
        self.rinj.unit = u.dimensionless_unscaled

        # Delta is of course unitless

        self.delta.unit = u.dimensionless_unscaled
        self.b.unit = u.dimensionless_unscaled

        # Piv has the same unit as energy (which is z)

        self.piv.unit = z_unit
        self.piv2.unit = z_unit

    def evaluate(self, x, y, z, lon0, lat0, rdiff0, rinj, delta, b, piv, piv2):

        lon, lat = x, y
        energy = z

        # energy in kev -> TeV.
        # NOTE: the use of piv2 is necessary to preserve dimensional correctness: the logarithm can only be taken
        # of a dimensionless quantity, so there must be a pivot there.

        e_energy_piv2 = 17. * \
            np.power(old_div(energy, piv2), 0.54 + 0.046 *
                     np.log10(old_div(energy, piv2)))
        e_piv_piv2 = 17. * \
            np.power(old_div(piv, piv2), 0.54 + 0.046 *
                     np.log10(old_div(piv, piv2)))

        rdiff_c = rdiff0 * np.power(old_div(e_energy_piv2, e_piv_piv2), old_div((delta - 1.), 2.)) * \
            np.sqrt(b * b / 8. / np.pi * 0.624 + 0.26 * np.power(1. + 0.0107 * e_piv_piv2, -1.5)) / \
            np.sqrt(b * b / 8. / np.pi * 0.624 + 0.26 *
                    np.power(1. + 0.0107 * e_energy_piv2, -1.5))

        rdiff_i = rdiff0 * rinj * \
            np.power(old_div(e_energy_piv2, e_piv_piv2), old_div(delta, 2.))

        rdiff = np.minimum(rdiff_c, rdiff_i)

        angsep = angular_distance_fast(lon, lat, lon0, lat0)

        pi = np.pi

        rdiffs, angseps = np.meshgrid(rdiff, angsep)

        return np.power(old_div(180.0, pi), 2) * 1.2154 / (pi * np.sqrt(pi) * rdiffs * (angseps + 0.06 * rdiffs)) * \
            np.exp(old_div(-np.power(angseps, 2), rdiffs ** 2))

    def get_boundaries(self):

        # Truncate the function at the max of rdiff allowed

        maximum_rdiff = self.rdiff0.max_value

        min_latitude = max(-90., self.lat0.value - maximum_rdiff)
        max_latitude = min(90., self.lat0.value + maximum_rdiff)

        max_abs_lat = max(np.absolute(min_latitude), np.absolute(max_latitude))

        if max_abs_lat > 89. or old_div(maximum_rdiff, np.cos(max_abs_lat * np.pi / 180.)) >= 180.:

            min_longitude = 0.
            max_longitude = 360.

        else:

            min_longitude = self.lon0.value - \
                old_div(maximum_rdiff, np.cos(max_abs_lat * np.pi / 180.))
            max_longitude = self.lon0.value + \
                old_div(maximum_rdiff, np.cos(max_abs_lat * np.pi / 180.))

            if min_longitude < 0.:

                min_longitude += 360.

            elif max_longitude > 360.:

                max_longitude -= 360.

        return (min_longitude, max_longitude), (min_latitude, max_latitude)

    def get_total_spatial_integral(self, z=None):
        """
        Returns the total integral (for 2D functions) or the integral over the spatial components (for 3D functions).
        needs to be implemented in subclasses.

        :return: an array of values of the integral (same dimension as z).
        """

        if isinstance(z, u.Quantity):
            z = z.value
        return np.ones_like(z)


class Continuous_injection_diffusion_legacy(Function3D, metaclass=FunctionMeta):
    r"""
        description :

            Positron and electrons diffusing away from the accelerator

        latex : $\left(\frac{180^\circ}{\pi}\right)^2 \frac{1.2154}{\sqrt{\pi^3} r_{\rm diff} ({\rm angsep} ({\rm x, y, lon_0, lat_0})+0.06 r_{\rm diff} )} \, {\rm exp}\left(-\frac{{\rm angsep}^2 ({\rm x, y, lon_0, lat_0})}{r_{\rm diff} ^2} \right)$

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

            rdiff0 :

                desc : Projected diffusion radius. The maximum allowed value is used to define the truncation radius.
                initial value : 1.0
                min : 0
                max : 20

            delta :

                desc : index for the diffusion coefficient
                initial value : 0.5
                min : 0.3
                max : 0.6
                fix : yes

            uratio :

                desc : ratio between u_cmb and u_B
                initial value : 0.5
                min : 0.01
                max : 100.
                fix : yes

            piv :

                desc : Pivot for the diffusion radius
                initial value : 2e10
                min : 0
                fix : yes

            piv2 :
                desc : Pivot for converting gamma energy to electron energy (always be 1 TeV)
                initial value : 1e9
                min : 0
                fix : yes

        """

    def _set_units(self, x_unit, y_unit, z_unit, w_unit):

        # lon0 and lat0 and rdiff have most probably all units of degrees. However,
        # let's set them up here just to save for the possibility of using the
        # formula with other units (although it is probably never going to happen)

        self.lon0.unit = x_unit
        self.lat0.unit = y_unit
        self.rdiff0.unit = x_unit

        # Delta is of course unitless

        self.delta.unit = u.dimensionless_unscaled
        self.uratio.unit = u.dimensionless_unscaled

        # Piv has the same unit as energy (which is z)

        self.piv.unit = z_unit
        self.piv2.unit = z_unit

    def evaluate(self, x, y, z, lon0, lat0, rdiff0, delta, uratio, piv, piv2):

        lon, lat = x, y
        energy = z

        # energy in kev -> TeV.
        # NOTE: the use of piv2 is necessary to preserve dimensional correctness: the logarithm can only be taken
        # of a dimensionless quantity, so there must be a pivot there.

        e_energy_piv2 = 17. * \
            np.power(old_div(energy, piv2), 0.54 + 0.046 *
                     np.log10(old_div(energy, piv2)))
        e_piv_piv2 = 17. * \
            np.power(old_div(piv, piv2), 0.54 + 0.046 *
                     np.log10(old_div(piv, piv2)))

        try:

            rdiff = rdiff0 * np.power(old_div(e_energy_piv2, e_piv_piv2), old_div((delta - 1.), 2.)) * \
                np.sqrt(1. + uratio * np.power(1. + 0.0107 * e_piv_piv2, -1.5)) / \
                np.sqrt(1. + uratio * np.power(1. +
                                               0.0107 * e_energy_piv2, -1.5))

        except ValueError:

            # This happens when using units, because astropy.units fails with the message:
            # "ValueError: Quantities and Units may only be raised to a scalar power"

            # Work around the problem with this loop, which is slow but using units is only for testing purposes or
            # single calls, so it shouldn't matter too much
            rdiff = np.array([(rdiff0 * np.power(old_div(e_energy_piv2, e_piv_piv2), x)).value for x in (delta - 1.) / 2. * np.sqrt(1. + uratio * np.power(1. + 0.0107 * e_piv_piv2, -1.5)) /
                              np.sqrt(1. + uratio * np.power(1. + 0.0107 * e_energy_piv2, -1.5))]) * rdiff0.unit

        angsep = angular_distance_fast(lon, lat, lon0, lat0)

        pi = np.pi

        rdiffs, angseps = np.meshgrid(rdiff, angsep)

        return np.power(old_div(180.0, pi), 2) * 1.2154 / (pi * np.sqrt(pi) * rdiffs * (angseps + 0.06 * rdiffs)) * \
            np.exp(old_div(-np.power(angseps, 2), rdiffs ** 2))

    def get_boundaries(self):

        # Truncate the function at the max of rdiff allowed

        maximum_rdiff = self.rdiff0.max_value

        min_latitude = max(-90., self.lat0.value - maximum_rdiff)
        max_latitude = min(90., self.lat0.value + maximum_rdiff)

        max_abs_lat = max(np.absolute(min_latitude), np.absolute(max_latitude))

        if max_abs_lat > 89. or old_div(maximum_rdiff, np.cos(max_abs_lat * np.pi / 180.)) >= 180.:

            min_longitude = 0.
            max_longitude = 360.

        else:

            min_longitude = self.lon0.value - \
                old_div(maximum_rdiff, np.cos(max_abs_lat * np.pi / 180.))
            max_longitude = self.lon0.value + \
                old_div(maximum_rdiff, np.cos(max_abs_lat * np.pi / 180.))

            if min_longitude < 0.:

                min_longitude += 360.

            elif max_longitude > 360.:

                max_longitude -= 360.

        return (min_longitude, max_longitude), (min_latitude, max_latitude)

    def get_total_spatial_integral(self, z=None):
        """
        Returns the total integral (for 2D functions) or the integral over the spatial components (for 3D functions).
        needs to be implemented in subclasses.

        :return: an array of values of the integral (same dimension as z).
        """

        if isinstance(z, u.Quantity):
            z = z.value
        return np.ones_like(z)


class GalPropTemplate_3D(Function3D):
    r"""
        description :

            Use a 3D template that has morphology and flux information.
            GalProp, DRAGON or a similar model in fits format would work. 
            Only parameter is a normalization factor. 

        latex : $ K $

        parameters :

            K :

                desc : normalization
                initial value : 1
                fix : yes

            hash :

                desc : hash of model map [needed for memoization]
                initial value : 1
                fix : yes         

    """

    __metaclass__ = FunctionMeta

    def _set_units(self, x_unit, y_unit, z_unit, w_unit):

        self.K.unit = (u.MeV * u.cm**2 * u.s * u.sr) ** (-1)

    def _setup(self):

        self._frame = ICRS()
        self._map = None
        self._fitsfile = None
        self._interpmap = None

    def set_frame(self, new_frame):
        """
        Set a new frame for the coordinates (the default is ICRS J2000)

        :param new_frame: a coordinate frame from astropy
        :return: (none)
        """
        assert isinstance(new_frame, BaseCoordinateFrame)

        self._frame = new_frame

    def load_file(self, fitsfile, phi1, phi2, theta1, theta2, galactic=False, ihdu=0):

        if fitsfile is None:
            raise RuntimeError(
                "Need to specify a fits file with a template map.")

        self._fitsfile = fitsfile
        p1, p2, t1, t2 = self.define_region(
            phi1, phi2, theta1, theta2, galactic)
        self.ramin = p1
        self.ramax = p2
        self.decmin = t1
        self.decmax = t2

        with fits.open(self._fitsfile) as f:

            self._delLon = f[ihdu].header['CDELT1']
            self._delLat = f[ihdu].header['CDELT2']
            self._delEn = f[ihdu].header['CDELT3']
            self._refLon = f[ihdu].header['CRVAL1']
            self._refLat = f[ihdu].header['CRVAL2']
            self._refEn = f[ihdu].header['CRVAL3']  # values in log10
            self._map = f[ihdu].data
            self._nl = f[ihdu].header['NAXIS1']  # longitude
            self._nb = f[ihdu].header['NAXIS2']  # latitude
            self._ne = f[ihdu].header['NAXIS3']  # energy

            # Create the function for the interpolation
            self._L = np.linspace(
                self._refLon, self._refLon+(self._nl-1)*self._delLon, self._nl)
            self._B = np.linspace(
                self._refLat, self._refLat+(self._nb-1)*self._delLat, self._nb)
            self._E = np.linspace(
                self._refEn, self._refEn+(self._ne-1)*self._delEn, self._ne)
            for i in range(len(self._E)):
                # Map units in Mev / cm^2 s sr, changing to 1 / MeV cm^2 s sr
                self._map[i] = self._map[i] / \
                    (np.power(10, self._E[i])*np.power(10, self._E[i]))
                self._map[i] = (np.fliplr(self._map[i]))
            self._F = RegularGridInterpolator(
                (self._E, self._B, self._L), self._map, bounds_error=False)

            h = hashlib.sha224()
            h.update(self._map)
            self.hash = int(h.hexdigest(), 16)

    def to_dict(self, minimal=False):

        data = super(Function3D, self).to_dict(minimal)

        if not minimal:

            data['extra_setup'] = {
                "_fitsfile": self._fitsfile,
                "_frame": self._frame,
                "ramin": self.ramin,
                "ramax": self.ramax,
                "decmin": self.decmin,
                "decmax": self.decmax
            }

        return data

    def which_model_file(self):
        return self._fitsfile

    def evaluate(self, x, y, z, K, hash):

        if self._map is None:

            self.load_file(self._fitsfile, self.ramin, self.ramax,
                           self.decmin, self.decmax, False, ihdu=0)

        # Interpolated values can be cached since we are fitting the constant K
        if self._interpmap is None:

            # We assume x and y are R.A. and Dec
            _coord = SkyCoord(ra=x, dec=y, frame=self._frame, unit="deg")
            b = _coord.transform_to('galactic').b.value
            l = _coord.transform_to('galactic').l.value
            lon = l
            lat = b

            # transform energy from keV to MeV. Galprop Model starts at 100 MeV
            energy = np.log10(z)-np.log10((u.MeV.to('keV')/u.keV).value)

            if lon.size != lat.size:
                raise AttributeError("Lon and Lat should be the same size")
            f = np.zeros([lon.size, energy.size])
            E0 = self._refEn
            Ef = self._refEn + (self._ne-1)*self._delEn

            # fix longitude
            shift = np.where(lon > 180.)
            lon[shift] = 180 - lon[shift]

            for i in range(energy.size):
                e = np.repeat(energy[i], len(lon))
                try:
                    f[:, i] = self._F(zip(e, lat, lon))
                except ValueError:
                    pass
            bad_idx = np.isnan(f)
            f[bad_idx] = 0
            bad_idx = np.isinf(f)
            f[bad_idx] = 0
            assert np.all(np.isfinite(f)), "some interpolated values are wrong"
            self._interpmap = f

        # (1000 is to change from MeV to KeV)
        A = np.multiply(K, self._interpmap/1000.)
        return A

    def define_region(self, a, b, c, d, galactic=False):
        if galactic:
            lmin = a
            lmax = b
            bmin = c
            bmax = d

            _coord = SkyCoord(l=[lmin, lmin, lmax, lmax], b=[
                              bmin, bmax, bmax, bmin], frame='galactic', unit='deg')

            ramin = min(_coord.transform_to('icrs').ra.value)
            ramax = max(_coord.transform_to('icrs').ra.value)
            decmin = min(_coord.transform_to('icrs').dec.value)
            decmax = max(_coord.transform_to('icrs').dec.value)

        else:
            ramin = a
            ramax = b
            decmin = c
            decmax = d
        return ramin, ramax, decmin, decmax

    def get_boundaries(self):
        min_longitude = self.ramin
        max_longitude = self.ramax
        min_latitude = self.decmin
        max_latitude = self.decmax
        return (min_longitude, max_longitude), (min_latitude, max_latitude)
