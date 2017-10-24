import numpy as np
import pdb

from scipy.interpolate import RegularGridInterpolator
from joblib import Parallel, delayed

from astropy.wcs import wcs
from astropy.coordinates import SkyCoord, ICRS, BaseCoordinateFrame
from astropy.io import fits
import astropy.units as u

from astromodels.functions.function import Function3D, FunctionMeta
from astromodels.utils.angular_distance import angular_distance



class Continuous_injection_diffusion_ellipse(Function3D):
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

        e_energy_piv2 = 17. * np.power(energy / piv2, 0.54 + 0.046 * np.log10(energy / piv2))
        e_piv_piv2 = 17. * np.power(piv / piv2, 0.54 + 0.046 * np.log10(piv / piv2))

        try:

            rdiff_a = rdiff0 * np.power(e_energy_piv2 / e_piv_piv2, (delta - 1.) / 2.) * \
                    np.sqrt(b * b / 8. / np.pi * 0.624 + 0.26 * np.power(1. + 0.0107 * e_piv_piv2, -1.5)) / \
                    np.sqrt(b * b / 8. / np.pi * 0.624 + 0.26 * np.power(1. + 0.0107 * e_energy_piv2, -1.5))

        except ValueError:

            # This happens when using units, because astropy.units fails with the message:
            # "ValueError: Quantities and Units may only be raised to a scalar power"

            # Work around the problem with this loop, which is slow but using units is only for testing purposes or
            # single calls, so it shouldn't matter too much
            rdiff_a = np.array( map(lambda x: (rdiff0 * np.power(e_energy_pivi2 / e_piv_piv2, x)).value,
                                  (delta - 1.) / 2. * np.sqrt(b * b / 8. / np.pi * 0.624 + 0.26 * np.power(1. + 0.0107 * e_piv_piv2, -1.5)) /
                                  np.sqrt(b * b / 8. / np.pi * 0.624 + 0.26 * np.power(1. + 0.0107 * e_energy_piv2, -1.5)))) * rdiff0.unit

        rdiff_b = rdiff_a * elongation

        pi = np.pi

        angsep = angular_distance(lon, lat, lon0, lat0)
        ang = np.arctan2(lat - lat0, (lon - lon0) * np.cos(lat0 * np.pi / 180.))

        theta = np.arctan2(np.sin(ang-incl*np.pi/180.)/elongation, np.cos(ang-incl*np.pi/180.))

        rdiffs_a, thetas = np.meshgrid(rdiff_a, theta)
        rdiffs_b, angseps = np.meshgrid(rdiff_b, angsep)

        rdiffs = np.sqrt(rdiffs_a ** 2 * np.cos(thetas) ** 2 + rdiffs_b ** 2 * np.sin(thetas) ** 2)


        results = np.power(180.0 / pi, 2) * 1.22 / (pi * np.sqrt(pi) * rdiffs_a * np.sqrt(elongation) * (angseps + 0.06 * rdiffs)) *  np.exp(-np.power(angseps, 2) / rdiffs ** 2)

        return results


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

        return (min_longitude, max_longitude), (min_latitude, max_latitude)

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

    __metaclass__ = FunctionMeta

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
        print z
        # energy in kev -> TeV.
        # NOTE: the use of piv2 is necessary to preserve dimensional correctness: the logarithm can only be taken
        # of a dimensionless quantity, so there must be a pivot there.

        e_energy_piv2 = 17. * np.power(energy / piv2, 0.54 + 0.046 * np.log10(energy / piv2))
        e_piv_piv2 = 17. * np.power(piv / piv2, 0.54 + 0.046 * np.log10(piv / piv2))

        try:

            rdiff_c = rdiff0 * np.power(e_energy_piv2 / e_piv_piv2, (delta - 1.) / 2.) * \
                    np.sqrt(b * b / 8. / np.pi * 0.624 + 0.26 * np.power(1. + 0.0107 * e_piv_piv2, -1.5)) / \
                    np.sqrt(b * b / 8. / np.pi * 0.624 + 0.26 * np.power(1. + 0.0107 * e_energy_piv2, -1.5))

            rdiff_i = rdiff0 * rinj * np.power(e_energy_piv2 / e_piv_piv2, delta / 2.)

        except ValueError:

            # This happens when using units, because astropy.units fails with the message:
            # "ValueError: Quantities and Units may only be raised to a scalar power"

            # Work around the problem with this loop, which is slow but using units is only for testing purposes or
            # single calls, so it shouldn't matter too much
            rdiff_c = np.array( map(lambda x: (rdiff0 * np.power(e_energy_pivi2 / e_piv_piv2, x)
                                   * np.sqrt(b * b / 8. / np.pi * 0.624 + 0.26 * np.power(1. + 0.0107 * e_piv_piv2, -1.5)) /
                                  np.sqrt(b * b / 8. / np.pi * 0.624 + 0.26 * np.power(1. + 0.0107 * e_energy_piv2, -1.5)), (delta - 1.) / 2.).value)) * rdiff0.unit

            rdiff_i = np.array( map(lambda x: (rdiff0 * rinj * np.power(e_energy_piv2 / e_piv_piv2, x), delta / 2.).value)) * rdiff0.unit
            raise ValueError("")

        rdiff = np.minimum(rdiff_c, rdiff_i)

        angsep = angular_distance(lon, lat, lon0, lat0)

        pi = np.pi

        rdiffs, angseps = np.meshgrid(rdiff, angsep)

        return np.power(180.0 / pi, 2) * 1.2154 / (pi * np.sqrt(pi) * rdiffs * (angseps + 0.06 * rdiffs)) * \
               np.exp(-np.power(angseps, 2) / rdiffs ** 2)


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

        return (min_longitude, max_longitude), (min_latitude, max_latitude)

class GalPropTemplate_3D(Function3D):
    r"""
        description :

            User input Spatial Template.

        latex : $ K $

        parameters :

            K :

                desc : normalization
                initial value : 1
                fix : yes

    """

    __metaclass__ = FunctionMeta

    def _set_units(self, x_unit, y_unit, z_unit, w_unit):

        self.K.unit = (u.MeV * u.cm**2 * u.s * u.sr) ** (-1)

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

    def load_file(self,fitsfile,ihdu=0):

        header = fits.getheader(fitsfile)
        self._w = wcs.WCS(header)
        self.ramin = 0.1
        self.ramax = 359.9
        self.decmin = -26.
        self.decmax = 66.
        with fits.open(fitsfile) as f:

            #self._refXpix = f[ihdu].header['CRPIX1']
            #self._refYpix = f[ihdu].header['CRPIX2']
            #self._refZpix = f[ihdu].header['CRPIX3']
            self._delLon = f[ihdu].header['CDELT1']
            self._delLat = f[ihdu].header['CDELT2']
            self._delEn = f[ihdu].header['CDELT3']
            self._refLon = f[ihdu].header['CRVAL1']
            self._refLat = f[ihdu].header['CRVAL2']
            self._refEn = f[ihdu].header['CRVAL3'] #in log10
            self._map = f[ihdu].data
            self._nl = f[ihdu].header['NAXIS1']#long
            self._nb = f[ihdu].header['NAXIS2']#lat
            self._ne = f[ihdu].header['NAXIS3']#energy

            #Create the function for the interpolation
            self._L = np.linspace(self._refLon,self._refLon+(self._nl-1)*self._delLon,self._nl)
            self._B = np.linspace(self._refLat,self._refLat+(self._nb-1)*self._delLat,self._nb)
            self._E = np.linspace(self._refEn,self._refEn+(self._ne-1)*self._delEn,self._ne)
            for i in xrange(len(self._E)):
                self._map[i] = self._map[i]/(np.power(10,self._E[i])*np.power(10,self._E[i])) # map is in Mev / cm^2 s sr, changing to 1 / MeV cm^2 s sr
            self._F = RegularGridInterpolator((self._E,self._B,self._L),self._map,bounds_error=False)

    def evaluate(self, x,y,z,K):

        # We assume x and y are R.A. and Dec
        _coord = SkyCoord(ra=x, dec=y, frame=self._frame, unit="deg")

        b = _coord.transform_to('galactic').b.value
        l = _coord.transform_to('galactic').l.value
        lon=l
        lat=b
        #transform energy from keV to MeV. Galprop Model starts at 100 MeV
        energy = np.log10(z * u.keV/ u.MeV)
        #print "Energies: ",np.power(10,energy)
        if lon.size != lat.size:
            raise AttributeError("Lon and Lat should be the same size")
        f=np.zeros([lon.size,energy.size])
        E0 = self._refEn
        Ef = self._refEn + (self._ne-1)*self._delEn

        #fix longitude
        for j in xrange(lon.size):
            if lon[j]>180.:
                lon[j]=180-lon[j]

        for i in xrange(energy.size):
            #print i
            if energy[i]<E0 or energy[i]>Ef:  #Maybe needed, it probably not necesary once the energy units are right?
                continue
            #r = Parallel(n_jobs=1)(delayed(self._F)((energy[i],lat[j],lon[j])) for j in xrange(lon.size))
            #f[:][i] = r

            for j in xrange(lon.size):
                try:
                    f[j,i] = self._F((energy[i],lat[j],lon[j]))#/(energy[i]*energy[i])
                except ValueError:
                    continue

        assert np.all(np.isfinite(f))
        A = np.multiply(K,f/1000.) #divide by 1000 since LiFF will mutliply it back (change from MeV to KeV)
        #print "Flux: ", A
        return A

    def define_region(self,a,b,c,d):
        if c < -26.:
            c = -26.
            print "Value cannot be lower than dec=-26 due to HAWC's FOV"
        if d > 66.:
            d = 66.
            print "Value cannot be lower than dec=66 due to HAWC's FOV"
        self.ramin = a
        self.ramax = b
        self.decmin = c
        self.decmax = d

    def get_boundaries(self):

        #if min_longitude < self.ramin:
        #    min_longitude = self.ramin
        #if max_longitude > self.ramax:
        #    max_longitude = self.ramax
        #if min_latitude < self.decmin:
        #    min_latitude = self.decmin
        #if max_latitude > self.decmax:
        #    max_latitude = self.decmax
        #maximum_rdiff = self.r

        #min_latitude = max(-90., self.lat0 - maximum_rdiff)
        #max_latitude = min(90., self.lat0 + maximum_rdiff)

        #max_abs_lat = max(np.absolute(min_latitude), np.absolute(max_latitude))

        #if max_abs_lat > 89. or maximum_rdiff / np.cos(max_abs_lat * np.pi / 180.) >= 180.:

           # min_longitude = 0.
           # max_longitude = 360.

        #else:

            #min_longitude = self.lon0 - maximum_rdiff / np.cos(max_abs_lat * np.pi / 180.)
            #max_longitude = self.lon0 + maximum_rdiff / np.cos(max_abs_lat * np.pi / 180.)

            #if min_longitude < 0.:

                #min_longitude += 360.

            #elif max_longitude > 360.:

                #max_longitude -= 360.
        #
        min_longitude = self.ramin
        max_longitude = self.ramax
        min_latitude = self.decmin
        max_latitude = self.decmax
        return (min_longitude, max_longitude), (min_latitude, max_latitude)
