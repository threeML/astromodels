import numpy as np
from scipy.interpolate import RegularGridInterpolator
import astropy.units as astropy_units

from astromodels.functions.function import Function1D, FunctionMeta
from astromodels.utils.data_files import _get_data_file_path


class DMFitFunction(Function1D):
    r"""
        description :

            Class that evaluates the spectrum for a DM particle of a given
            mass, channel, cross section, and J-factor. Based on standard
            Fermi Science Tools function DMFitFunction. Note input table only
            calculated spectra up to m_DM of 10 TeV

            The parameterization is given by

            F(x) = 1 / (8 * pi) * (1/mass^2) * sigmav * J * dN/dE(E,mass,i)

        latex : $$

        parameters :

            mass :
                desc : DM mass (GeV)
                initial value : 10
                fix : yes

            channel :
                desc : DM annihilation channel
                initial value : 4
                fix : yes

            sigmav :
                desc : DM annihilation cross section (cm^3/s)
                initial value : 1.e-26

            J :
                desc : Target total J-factor (GeV^2 cm^-5)
                initial value : 1.e20
                fix : yes
        """

    __metaclass__ = FunctionMeta

    def _setup(self):

        tablepath = _get_data_file_path("dark_matter/gammamc_dif.dat")

        self._data = np.loadtxt(tablepath)

        """
            Mapping between the channel codes and the rows in the gammamc file

            1 : 8, # ee
            2 : 6, # mumu
            3 : 3, # tautau
            4 : 1, # bb
            5 : 2, # tt
            6 : 7, # gg
            7 : 4, # ww
            8 : 5, # zz
            9 : 0, # cc
            10 : 10, # uu
            11 : 11, # dd
            12 : 9, # ss
        """

        channel_index_mapping = {
            1: 8,  # ee
            2: 6,  # mumu
            3: 3,  # tautau
            4: 1,  # bb
            5: 2,  # tt
            6: 7,  # gg
            7: 4,  # ww
            8: 5,  # zz
            9: 0,  # cc
            10: 10,  # uu
            11: 11,  # dd
            12: 9,  # ss
        }

        # Number of decades in x = log10(E/M)
        ndec = 10.0
        xedge = np.linspace(0, 1.0, 251)
        self._x = 0.5 * (xedge[1:] + xedge[:-1]) * ndec - ndec

        ichan = channel_index_mapping[int(self.channel.value)]

        # These are the mass points
        self._mass = np.array([2.0, 4.0, 6.0, 8.0, 10.0,
                               25.0, 50.0, 80.3, 91.2, 100.0,
                               150.0, 176.0, 200.0, 250.0, 350.0, 500.0, 750.0,
                               1000.0, 1500.0, 2000.0, 3000.0, 5000.0, 7000.0, 1E4])
        self._dn = self._data.reshape((12, 24, 250))

        self._dn_interp = RegularGridInterpolator([self._mass, self._x],
                                                  self._dn[ichan, :, :],
                                                  bounds_error=False,
                                                  fill_value=None)

        if self.mass.value > 10000:

            print "Warning: DMFitFunction only appropriate for masses <= 10 TeV"
            print "To model DM from 2 GeV < mass < 1 PeV use DMSpectra"

    def _set_units(self, x_unit, y_unit):

        # Usually a model should not assume fixed units for energy or anything else. However,
        # in this case this model is so specialistic that we can assume GeV

        self.mass.unit = astropy_units.GeV

        self.channel.unit = astropy_units.dimensionless_unscaled

        self.sigmav.unit = astropy_units.cm ** 3 / astropy_units.s

        self.J.unit = astropy_units.GeV ** 2 / astropy_units.cm ** 5

    def print_channel_mapping(self):

        channel_mapping = {
            1: 'ee',
            2: 'mumu',
            3: 'tautau',
            4: 'bb',
            5: 'tt',
            6: 'gg',
            7: 'ww',
            8: 'zz',
            9: 'cc',
            10: 'uu',
            11: 'dd',
            12: 'ss',
        }

        print channel_mapping

        return channel_mapping

    # noinspection PyPep8Naming
    def evaluate(self, x, mass, channel, sigmav, J):

        if isinstance(x, astropy_units.Quantity):

            # We need to convert to GeV
            xx = x.to(astropy_units.GeV)

        else:

            # We can assume that the input is in keV

            keVtoGeV = 1e-6

            xx = np.multiply(x, keVtoGeV)  # xm expects gamma ray energies in MeV

        xm = np.log10(np.divide(xx, mass))
        phip = 1. / (8. * np.pi) * np.power(mass, -2) * (sigmav * J)  # units of this should be 1 / cm**2 / s
        dn = self._dn_interp((mass, xm))
        dn[xm > 0] = 0

        return np.multiply(phip, np.divide(dn, x))


class DMSpectra(Function1D):
    r"""
        description :

            Class that evaluates the spectrum for a DM particle of a given
            mass, channel, cross section, and J-factor. Combines Pythia-based tables
            from both Fermi (2 GeV < m_DM < 10 TeV) and HAWC (10 TeV < m_dm < 1 PeV)

            The parameterization is given by

            F(x) = 1 / (8 * pi) * (1/mass^2) * sigmav * J * dN/dE(E,mass,i)

            Note that this class assumes that mass and J-factor are provided
            in units of GeV and GeV^2 cm^-5

        latex : $$

        parameters :

            mass :
                desc : DM mass (GeV)
                initial value : 10
                fix : yes

            channel :
                desc : DM annihilation channel
                initial value : 4
                fix : yes

            sigmav :
                desc : DM annihilation cross section (cm^3/s)
                initial value : 1.e-26

            J :
                desc : Target total J-factor (GeV^2 cm^-5)
                initial value : 1.e20
                fix : yes
        """

    __metaclass__ = FunctionMeta

    def _setup(self):

        # Get and open the two data files

        tablepath_h = _get_data_file_path("dark_matter/dmSpecTab.npy")
        self._data_h = np.load(tablepath_h)

        tablepath_f = _get_data_file_path("dark_matter/gammamc_dif.dat")
        self._data_f = np.loadtxt(tablepath_f)

        """
            Mapping between the channel codes and the rows in the gammamc file
            dmSpecTab.npy created to match this mapping too

            1 : 8, # ee
            2 : 6, # mumu
            3 : 3, # tautau
            4 : 1, # bb
            5 : 2, # tt
            6 : 7, # gg
            7 : 4, # ww
            8 : 5, # zz
            9 : 0, # cc
            10 : 10, # uu
            11 : 11, # dd
            12 : 9, # ss
            """

        channel_index_mapping = {
            1: 8,  # ee
            2: 6,  # mumu
            3: 3,  # tautau
            4: 1,  # bb
            5: 2,  # tt
            6: 7,  # gg
            7: 4,  # ww
            8: 5,  # zz
            9: 0,  # cc
            10: 10,  # uu
            11: 11,  # dd
            12: 9,  # ss
        }

        # Number of decades in x = log10(E/M)
        ndec = 10.0
        xedge = np.linspace(0, 1.0, 251)
        self._x = 0.5 * (xedge[1:] + xedge[:-1]) * ndec - ndec

        ichan = channel_index_mapping[int(self.channel.value)]

        # These are the mass points in GeV
        self._mass_h = np.array([50., 61.2, 74.91, 91.69, 112.22, 137.36, 168.12, 205.78, 251.87, 308.29,
                                 377.34, 461.86, 565.31, 691.93, 846.91, 1036.6, 1268.78, 1552.97, 1900.82,
                                 2326.57, 2847.69, 3485.53, 4266.23, 5221.81, 6391.41, 7823.0, 9575.23,
                                 11719.94, 14345.03, 17558.1, 21490.85, 26304.48, 32196.3, 39407.79, 48234.54,
                                 59038.36, 72262.07, 88447.7, 108258.66, 132506.99, 162186.57, 198513.95,
                                 242978.11, 297401.58, 364015.09, 445549.04, 545345.37, 667494.6, 817003.43, 1000000.])

        # These are the mass points in GeV
        self._mass_f = np.array([2.0, 4.0, 6.0, 8.0, 10.0,
                                 25.0, 50.0, 80.3, 91.2, 100.0,
                                 150.0, 176.0, 200.0, 250.0, 350.0, 500.0, 750.0,
                                 1000.0, 1500.0, 2000.0, 3000.0, 5000.0, 7000.0, 1E4])

        self._mass = np.append(self._mass_f, self._mass_h[27:])

        self._dn_f = self._data_f.reshape((12, 24, 250))

        # Is this really used?
        self._dn_h = self._data_h

        self._dn = np.zeros((12, len(self._mass), 250))
        self._dn[:, 0:24, :] = self._dn_f
        self._dn[:, 24:, :] = self._dn_h[:, 27:, :]

        self._dn_interp = RegularGridInterpolator([self._mass, self._x],
                                                  self._dn[ichan, :, :],
                                                  bounds_error=False,
                                                  fill_value=None)

        if self.channel.value in [1, 6, 7] and self.mass.value > 10000.:
            print "ERROR: currently spectra for selected channel and mass not implemented."
            print "Spectra for channels ['ee','gg','WW'] currently not available for mass > 10 TeV"

    def _set_units(self, x_unit, y_unit):

        self.mass.unit = astropy_units.GeV
        self.channel.unit = astropy_units.dimensionless_unscaled
        self.sigmav.unit = astropy_units.cm ** 3 / astropy_units.s
        self.J.unit = astropy_units.GeV ** 2 / astropy_units.cm ** 5

    def print_channel_mapping(self):
        channel_mapping = {
            1: 'ee',
            2: 'mumu',
            3: 'tautau',
            4: 'bb',
            5: 'tt',
            6: 'gg',
            7: 'ww',
            8: 'zz',
            9: 'cc',
            10: 'uu',
            11: 'dd',
            12: 'ss',
        }

        print channel_mapping

        return channel_mapping

    # noinspection PyPep8Naming
    def evaluate(self, x, mass, channel, sigmav, J):

        if isinstance(x, astropy_units.Quantity):

            # We need to convert to GeV
            xx = x.to(astropy_units.MeV)

        else:

            # We can assume that the input is in keV

            keVtoGeV = 1E-6

            xx = np.multiply(x, keVtoGeV)  # xm expects gamma ray energies in MeV

        xm = np.log10(np.divide(xx, mass))

        phip = 1. / (8. * np.pi) * np.power(mass, -2) * (sigmav * J)  # units of this should be 1 / cm**2
        dn = self._dn_interp((mass, xm))  # note this is unitless (dx = d(xm))
        dn[xm > 0] = 0

        return np.multiply(phip, np.divide(dn, x))
