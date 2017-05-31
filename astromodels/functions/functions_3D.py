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
        self.uratio.unit = u.dimensionless_unscaled

        # Piv has the same unit as energy (which is z)

        self.piv.unit = z_unit
        self.piv2.unit = z_unit

    def evaluate(self, x, y, z, lon0, lat0, rdiff0, delta, uratio, piv, piv2):

        print x
        print y
        print z.size
        lon, lat = x, y
        energy = z

        # energy in kev -> TeV.
        # NOTE: the use of piv2 is necessary to preserve dimensional correctness: the logarithm can only be taken
        # of a dimensionless quantity, so there must be a pivot there.

        e_energy_piv2 = 17. * np.power(energy / piv2, 0.54 + 0.046 * np.log10(energy / piv2))
        e_piv_piv2 = 17. * np.power(piv / piv2, 0.54 + 0.046 * np.log10(piv / piv2))

        try:

            rdiff = rdiff0 * np.power(e_energy_piv2 / e_piv_piv2, (delta - 1.) / 2.) * \
                    np.sqrt(1. + uratio * np.power(1. + 0.0107 * e_piv_piv2, -1.5)) / \
                    np.sqrt(1. + uratio * np.power(1. + 0.0107 * e_energy_piv2, -1.5))

        except ValueError:

            # This happens when using units, because astropy.units fails with the message:
            # "ValueError: Quantities and Units may only be raised to a scalar power"

            # Work around the problem with this loop, which is slow but using units is only for testing purposes or
            # single calls, so it shouldn't matter too much
            rdiff = np.array( map(lambda x: (rdiff0 * np.power(e_energy_pivi2 / e_piv_piv2, x)).value,
                                  (delta - 1.) / 2. * np.sqrt(1. + uratio * np.power(1. + 0.0107 * e_piv_piv2, -1.5)) /
                                  np.sqrt(1. + uratio * np.power(1. + 0.0107 * e_energy_piv2, -1.5)))) * rdiff0.unit

        angsep = angular_distance(lon, lat, lon0, lat0)

        pi = np.pi

        rdiffs, angseps = np.meshgrid(rdiff, angsep)

        A = np.power(180.0 / pi, 2) * 1.2154 / (pi * np.sqrt(pi) * rdiffs * (angseps + 0.06 * rdiffs)) * \
               np.exp(-np.power(angseps, 2) / rdiffs ** 2)
        print A
        return A


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
            self._F = RegularGridInterpolator((self._E,self._B,self._L),self._map,bounds_error=False) 

    #deprecated, delete after testing the other one
    #def _interpolate_method(self,i,j,k,l,b,e):
        #energy pixels
        #k1 = int(np.floor(k))
        #k2 = int(np.ceil(k))

        #lon pixels
        #i1 = int(np.floor(i))
        #i2 = int(np.ceil(i))

        #lat pixels
        #j1 = int(np.floor(j))
        #j2 = int(np.ceil(j))

        #if in the edge of the map, use the value of the closest model pixel
        #if k2 == self._ne:
            #w1 = 1
            #w2 = 0
        #else:
            #Linear interoplation for energy
            #E1 = self._refEn + k1*self._delEn
            #E2 = self._refEn + k2*self._delEn
            #w1 = (e - E1)/self._delEn
            #w2 = (E2 - e)/self._delEn
        #f1 = self._map[k1][j1][i1]
        #f2 = self._map[k2][j1][i1]

        #if j2 == self._nb or i2 == self._nl:
        #    f1 = self._map[k1][j1][i1]
        #    if k2 == self._ne:
        #        f2 = 0.
        #    else:
        #        f2 = self._map[k2][j1][i1]
        #elif j1 == -1 or i1 == -1:
        #    f1 = self._map[k1][j2][i2]
        #    if k2 == self._ne:
        #        f2 = 0
        #    else:
        #        f2 = self._map[k2][j2][i2]
        #else: #do full bilinear interpolation
        #    y1 = self._refLat + self._delLat*j1
        #    x1 = self._refLon + self._delLon*i1
        #    y2 = self._refLat + self._delLat*j2
        #    x2 = self._refLon + self._delLon*i2
        #    ## Bilinear interpolation
        #    ##Getting flux from k1
        #    Q11 = self._map[k1][j1][i1]
        #    Q21 = self._map[k1][j1][i2]
        #    Q12 = self._map[k1][j2][i1]
        #    Q22 = self._map[k1][j2][i2]
        #    f1 = (Q11*(x2-l)*(y2-b) + Q21*(l-x1)*(y2-b) + Q12*(x2-l)*(b-y1) + Q22*(l-x2)*(b-y2))/(self._delLon*self._delLat)

        #    ##Getting flux from k2
        #    if k2 == self._ne:
        #        f2 = 0.
        #    else:
        #        Q11 = self._map[k2][j1][i1]
        #        Q21 = self._map[k2][j1][i2]
        #        Q12 = self._map[k2][j2][i1]
        #        Q22 = self._map[k2][j2][i2]
        #        f2 = (Q11*(x2-l)*(y2-b) + Q21*(l-x1)*(y2-b) + Q12*(x2-l)*(b-y1) + Q22*(l-x2)*(b-y2))/(self._delLon*self._delLat)
        #return w1*f1 + w2*f2

    def evaluate(self, x,y,z,K):

        # We assume x and y are R.A. and Dec
        print np.min(x)
        print x
        print y
        print x.size
        print z.size
        _coord = SkyCoord(ra=x, dec=y, frame=self._frame, unit="deg")

        b = _coord.transform_to('galactic').b.value
        l = _coord.transform_to('galactic').l.value
        lon=l
        lat=b
        #transform energy from keV to MeV. Galprop Model starts at 100 MeV
        #print np.log10(1.* u.keV/u.eV)
        energy = np.log10(z /1000.)#* u.eV/ u.MeV)
        print energy
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
                #f[j,i] = r[j]
                #il,ib,ie = self._w.all_world2pix(lon[j],lat[j],energy[i],1)
                #if il > self._nl+1:
                    #continue#il = il - self._nl
                #if ib > self._nb:
                    #ib = ib - self._nb
                #if ie > self._ne or energy[i]>Ef:  #Maybe needed, it probably not necesary once the energy units are right?
            #    if energy[i]<E0 or energy[i]>Ef:  #Maybe needed, it probably not necesary once the energy units are right?
                    #f[j,i] = 0.
                    #ie=ie-1.
                    #print energy[i]
                    #print "Looking at energies higher than 100 TeV"
             #       continue#f[j,i] = self._interpolate_method(il,ib,self._ne-1,lon[j],lat[j],energy[i])

                #else:
                    #f[j,i] = self._interpolate_method(il,ib,ie,lon[j],lat[j],energy[i])

                try:
                    f[j,i] = self._F((energy[i],lat[j],lon[j]))
                except ValueError:
                    continue

        pdb.set_trace()
        assert np.all(np.isfinite(f))
        A = np.multiply(K,f)
        print A
        return A

    def define_region(self,r,b0,l0):
        self.r = r
        self.lat0 = b0
        self.lon0 = l0
    
    def get_boundaries(self):
        maximum_rdiff = self.r

        min_latitude = max(-90., self.lat0 - maximum_rdiff)
        max_latitude = min(90., self.lat0 + maximum_rdiff)

        max_abs_lat = max(np.absolute(min_latitude), np.absolute(max_latitude))

        if max_abs_lat > 89. or maximum_rdiff / np.cos(max_abs_lat * np.pi / 180.) >= 180.:

            min_longitude = 0.
            max_longitude = 360.

        else:

            min_longitude = self.lon0 - maximum_rdiff / np.cos(max_abs_lat * np.pi / 180.)
            max_longitude = self.lon0 + maximum_rdiff / np.cos(max_abs_lat * np.pi / 180.)

            if min_longitude < 0.:

                min_longitude += 360.

            elif max_longitude > 360.:

                max_longitude -= 360.

        return (min_longitude, max_longitude), (min_latitude, max_latitude)
