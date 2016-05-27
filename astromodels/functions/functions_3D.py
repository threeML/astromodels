import astropy.units as u

from astromodels.functions.function import Function3D, FunctionMeta

import numpy as np

from astromodels.utils.angular_distance import angular_distance


class Continuous_injection_diffusion(Function3D):
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

            piv :

                desc : Pivot
                initial value : 2e10
                min : 0
                fix : yes

            piv2 :
                desc : Pivot for the diffusion radius
                initial value : 1e9
                min : 0
                fix : yes

        """

    __metaclass__ = FunctionMeta

    def _set_units(self, x_unit, y_unit, z_unit, w_unit):

        # lon0 and lat0 and rdiff have most probably all units of degrees. However,
        # let's set them up here just to save for the possibility of using the
        # formula with other units (although it is probably never going to happen)

        self.lon0.unit = x_unit
        self.lat0.unit = y_unit
        self.rdiff0.unit = x_unit

        # Delta is of course unitless

        self.delta.unit = u.dimensionless_unscaled

        # Piv has the same unit as energy (which is z)

        self.piv.unit = z_unit
        self.piv2.unit = z_unit

    def evaluate(self, x, y, z, lon0, lat0, rdiff0, delta, piv, piv2):

        lon, lat = x, y
        energy = z

        # energy in kev -> TeV.
        # NOTE: the use of piv2 is necessary to preserve dimensional correctness: the logarithm can only be taken
        # of a dimensionless quantity, so there must be a pivot there.
        try:

            rdiff = rdiff0 * np.power(energy / piv, (delta - 1.) / 2. * (0.54 + 0.046 * np.log10(energy / piv2)))

        except ValueError:

            # This happens when using units, because astropy.units fails with the message:
            # "ValueError: Quantities and Units may only be raised to a scalar power"

            # Work around the problem with this loop, which is slow but using units is only for testing purposes or
            # single calls, so it shouldn't matter too much
            rdiff = np.array( map(lambda x: (rdiff0 * np.power(energy/ piv, x)).value,
                                  (delta - 1.) / 2. * (0.54 + 0.046 * np.log10(energy / piv2)))) * rdiff0.unit

        angsep = angular_distance(lon, lat, lon0, lat0)

        pi = np.pi

        return np.power(180.0 / pi, 2) * 1.2154 / (pi * np.sqrt(pi) * rdiff * (angsep + 0.06 * rdiff)) * \
               np.exp(-np.power(angsep, 2) / rdiff ** 2)


    def get_boundaries(self):

        # Truncate the function at the max of rdiff allowed

        maximum_rdiff = self.rdiff0.max_value

        min_latitude = max(-90., self.lat0.value - maximum_rdiff)
        max_latitude = min(90., self.lat0.value + maximum_rdiff)

        max_abs_lat = max(np.absolute(min_latitude), np.absolute(max_latitude))

        if max_abs_lat > 89. or maximum_rdiff / np.cos(max_abs_lat * np.pi / 180.) >= 180.:

            min_longitude = 0.
            max_longitude = 360.

        else:

            min_longitude = self.lon0.value - maximum_rdiff / np.cos(max_abs_lat * np.pi / 180.)
            max_longitude = self.lon0.value + maximum_rdiff / np.cos(max_abs_lat * np.pi / 180.)

            if min_longitude < 0.:

                min_longitude += 360.

            elif max_longitude > 360.:

                max_longitude -= 360.

        return (min_longitude, max_longitude), (min_latitude, max_latitude), (0, np.inf)