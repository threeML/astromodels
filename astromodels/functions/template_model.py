from astromodels.functions.function import Function1D, FunctionMeta
from astromodels.utils.configuration import get_user_data_path

import collections
from astromodels.parameter import Parameter
import numpy as np
import pandas as pd
from pandas import HDFStore
import scipy.interpolate
import os
import re
import warnings


class IncompleteGrid(RuntimeError):
    pass


class ValuesNotInGrid(ValueError):
    pass


class MissingDataFile(RuntimeError):
    pass


class TemplateModelFactory(object):

    def __init__(self, name, description, energies, names_of_parameters,
                 interpolation_degree=1, flux_log=True,
                 spline_smoothing_factor=1,
                 nuFnu=True):

        # Store model name

        # Enforce that it does not contain spaces nor strange characters

        name = str(name)

        if re.match("[a-zA-Z_][a-zA-Z0-9_]*", name) is None:

            raise RuntimeError("The provided name '%s' is not a valid name. You cannot use spaces, "
                               "or special characters")

        self._name = name

        self._description = str(description)

        # Store energy grid

        self._energies = np.array(energies)

        # Enforce that they are ordered
        self._energies.sort()

        # We create a dictionary which will contain the grid for each parameter

        self._parameters_grids = collections.OrderedDict()

        for parameter_name in names_of_parameters:

            self._parameters_grids[parameter_name] = None

        self._data_frame = None
        self._multi_index = None
        self._interpolators = None

        self._interpolation_degree = interpolation_degree

        self._flux_log = bool(flux_log)

        self._spline_smoothing_factor = int(spline_smoothing_factor)

        # This indicate that the data which will be given contains the nuFnu spectrum
        self._nuFnu = bool(nuFnu)

    def define_parameter_grid(self, parameter_name, grid):

        assert parameter_name in self._parameters_grids, "Parameter %s is not part of this model" % parameter_name

        grid_ = np.array(grid)

        assert grid_.shape[0] > 1, "A grid for a parameter must contain at least two elements"

        # Assert that elements are unique

        assert np.all(np.unique(grid_) == grid_), "Non-unique elements in grid for parameter %s" % parameter_name

        self._parameters_grids[parameter_name] = grid_

    def add_interpolation_data(self, parameters_values, differential_fluxes):

        # Verify that the grid has been defined for all parameters

        for grid in self._parameters_grids.values():

            if grid is None:

                raise IncompleteGrid("You need to define a grid for all parameters, by using the "
                                     "define_parameter_grid method.")

        if self._data_frame is None:

            # This is the first data set, create the data frame

            # Create the multi-index

            self._multi_index = pd.MultiIndex.from_product(self._parameters_grids.values(),
                                                           names=self._parameters_grids.keys())

            # Pre-fill the data matrix with nans, so we will know if some elements have not been filled

            self._data_frame = pd.DataFrame(index=self._multi_index, columns=self._energies)

        # Make sure we are dealing with arrays (list will be transformed)

        parameters_values = np.array(parameters_values)
        differential_fluxes = np.array(differential_fluxes)

        n_parameters = parameters_values.shape[0]

        assert n_parameters == len(self._parameters_grids), "Wrong number of parameters"

        assert self._energies.shape[0] == differential_fluxes.shape[0], "Differential fluxes and energies must have " \
                                                                        "the same number of elements"

        # Now set the corresponding values in the data frame

        # Take the log of the flux, if requested

        if self._flux_log:

            values = np.log10(differential_fluxes)

        else:

            values = differential_fluxes

        if self._flux_log:

            msg = 'log10 of '

        else:

            msg = ''

        assert np.all(np.isfinite(values)), "One or more inf. or nan in the %svalues" % msg

        # Now set the values in the data frame

        try:

            self._data_frame.loc[tuple(parameters_values)] = pd.to_numeric(values)

        except KeyError:

            raise ValuesNotInGrid("The provided parameter values (%s) are not in the defined grid" % parameters_values)

    @staticmethod
    def _clean_cols_for_hdf(data):

        types = data.apply(lambda x: pd.lib.infer_dtype(x.values))

        for col in types.index:

            data[col] = pd.to_numeric(data[col])

        return data

    def save_data(self, overwrite=False):

        # First make sure that the whole data matrix has been filled

        assert not self._data_frame.isnull().values.any(), "You have NaNs in the data matrix. Usually this means " \
                                                           "that you didn't fill it up completely, or that some of " \
                                                           "your data contains nans. Cannot save the file."

        # Get the data directory

        data_dir_path = get_user_data_path()

        # Sanitize the data file

        filename_sanitized = os.path.abspath(os.path.join(data_dir_path, '%s.h5' % self._name))

        # Check that it does not exists
        if os.path.exists(filename_sanitized):

            if overwrite:

                try:

                    os.remove(filename_sanitized)

                except:

                    raise IOError("The file %s already exists and cannot be removed (maybe you do not have "
                                  "permissions to do so?). " % filename_sanitized)

            else:

                raise IOError("The file %s already exists! You cannot call two different "
                              "template models with the same name" % filename_sanitized)

        # Open the HDF5 file and write objects

        with HDFStore(filename_sanitized) as store:

            # The _clean_cols_for_hdf is needed because for some reasons the format of some columns
            # is not accepted by .to_hdf otherwise

            self._clean_cols_for_hdf(self._data_frame).to_hdf(store, 'data_frame')

            store.get_storer('data_frame').attrs.metadata = {'description': self._description,
                                                             'name': self._name,
                                                             'interpolation_degree': int(self._interpolation_degree),
                                                             'flux_log': self._flux_log,
                                                             'spline_smoothing_factor': self._spline_smoothing_factor,
                                                             'nuFnu': self._nuFnu}

            for i, parameter_name in enumerate(self._parameters_grids.keys()):

                store['p_%i_%s' % (i, parameter_name)] = pd.Series(self._parameters_grids[parameter_name])

            store['energies'] = pd.Series(self._energies)

# This adds a method to a class at runtime

def add_method(self, method, name=None):

    if name is None:

        name = method.func_name

    setattr(self.__class__, name, method)


class RectBivariateSplineWrapper(object):
    """
    Wrapper around RectBivariateSpline, which supplies a __call__ method which accept the same
    syntax as the other interpolation methods

    """

    def __init__(self, *args, **kwargs):

        # We can use interp2, which features spline interpolation instead of linear interpolation

        self._interpolator = scipy.interpolate.RectBivariateSpline(*args, **kwargs)

    def __call__(self, x):

        res = self._interpolator(*x)

        return res[0][0]


class TemplateModel(Function1D):

    r"""
        description :

            A template model

        latex : $n.a.$

        parameters :

            K :

                desc : Normalization (freeze this to 1 if the template provides the normalization by itself)
                initial value : 1.0

            scale :

                desc : Scale for the independent variable. The templates are handled as if they contains the fluxes
                       at x / scale. This is useful for example when the template describe energies in the rest frame,
                       at which point the scale describe the transformation between rest frame energy and observer frame
                       energy. Fix this to 1 to neutralize its effect.

                initial value : 1.0
                min : 1e-5

        """

    __metaclass__ = FunctionMeta

    def _custom_init_(self, model_name):

        # Get the data directory

        data_dir_path = get_user_data_path()

        # Sanitize the data file

        filename_sanitized = os.path.abspath(os.path.join(data_dir_path, '%s.h5' % model_name))

        if not os.path.exists(filename_sanitized):

            raise MissingDataFile("The data file %s does not exists. Did you use the "
                                  "TemplateFactory?" % (filename_sanitized))

        # Open the template definition and read from it

        self._data_file = filename_sanitized

        with HDFStore(filename_sanitized) as store:

            self._data_frame = store['data_frame']

            self._parameters_grids = collections.OrderedDict()

            processed_parameters = 0

            for key in store.keys():

                match = re.search('p_([0-9]+)_(.+)', key)

                if match is None:

                    continue

                else:

                    tokens = match.groups()

                    this_parameter_number = int(tokens[0])
                    this_parameter_name = str(tokens[1])

                    assert this_parameter_number == processed_parameters, "Parameters out of order!"

                    self._parameters_grids[this_parameter_name] = store[key]

                    processed_parameters += 1

            self._energies = store['energies']

            # Now get the metadata

            metadata = store.get_storer('data_frame').attrs.metadata

            description = metadata['description']
            name = metadata['name']
            self._flux_log = metadata['flux_log']

            self._interpolation_degree = metadata['interpolation_degree']

            self._spline_smoothing_factor = metadata['spline_smoothing_factor']

            self._nuFnu = metadata['nuFnu']

        # Make the dictionary of parameters

        function_definition = collections.OrderedDict()

        function_definition['description'] = description

        function_definition['latex'] = 'n.a.'

        # Now build the parameters according to the content of the parameter grid

        parameters = collections.OrderedDict()

        parameters['K'] = Parameter('K', 1.0)
        parameters['scale'] = Parameter('scale', 1.0)

        for parameter_name in self._parameters_grids.keys():

            grid = self._parameters_grids[parameter_name]

            parameters[parameter_name] = Parameter(parameter_name, grid.median(),
                                                   min_value=grid.min(),
                                                   max_value=grid.max())

        super(TemplateModel, self).__init__(name, function_definition, parameters)

        self._prepare_interpolators()

        # Now susbsitute the evaluate function with a version with all the required parameters

        # Get the parameters' names (except for K and scale)
        par_names_no_K_no_scale = parameters.keys()[2:]

        function_code = 'def new_evaluate(self, x, %s): ' \
                        'return K * self._interpolate(x, scale, [%s])' % (",".join(parameters.keys()),
                                                                          ",".join(par_names_no_K_no_scale))

        exec(function_code)

        add_method(self, new_evaluate,'_evaluate')

        self.evaluate = self._evaluate

    def _prepare_interpolators(self):

        # Figure out the shape of the data matrices
        data_shape = map(lambda x: x.shape[0], self._parameters_grids.values())

        self._interpolators = []

        for energy in self._energies:

            # Make interpolator for this energy
            this_data = np.array(self._data_frame[energy].values.reshape(*data_shape),dtype=float)

            if len(self._parameters_grids.values()) == 2:

                x, y = self._parameters_grids.values()

                # Make sure that the requested polynomial degree is less than the number of data sets in
                # both directions

                msg = "You cannot use an interpolation degree of %s if you don't provide at least %s points " \
                      "in the %s direction. Increase the number of templates or decrease the interpolation " \
                      "degree."

                if len(x) <= self._interpolation_degree:

                    raise RuntimeError(msg % (self._interpolation_degree, self._interpolation_degree+1, 'x'))

                if len(y) <= self._interpolation_degree:

                    raise RuntimeError(msg % (self._interpolation_degree, self._interpolation_degree + 1, 'y'))


                this_interpolator = RectBivariateSplineWrapper(x, y, this_data,
                                                               kx=self._interpolation_degree,
                                                               ky=self._interpolation_degree,
                                                               s=self._spline_smoothing_factor)

            else:

                # In more than 2d we can only use linear interpolation

                this_interpolator = scipy.interpolate.RegularGridInterpolator(self._parameters_grids.values(),
                                                                              this_data)

            self._interpolators.append(this_interpolator)

    def _set_units(self, x_unit, y_unit):

        pass

    # This function will be substituted during construction by another version with
    # all the parameters of this template

    def evaluate(self, x, K, scale):

        raise NotImplementedError("Should not get here!")

    def _interpolate(self, energies, scale, parameters_values):

        log_energies = np.log10(energies)

        e_tilde = self._energies * scale

        # Gather all interpolations for these parameters' values at all defined energies

        interpolations = np.array(map(lambda i:self._interpolators[i](parameters_values),
                                      range(self._energies.shape[0])))

        # Now interpolate the interpolations to get the flux at the requested energies

        if self._flux_log:

            # NOTE: the variable "interpolations" contains already the log10 of the values,
            # if _flux_log is True

            interpolator = scipy.interpolate.InterpolatedUnivariateSpline(np.log10(e_tilde),
                                                                          interpolations,
                                                                          k=self._interpolation_degree,
                                                                          ext=0)

            values = np.power(10, interpolator(log_energies))

        else:

            interpolator = scipy.interpolate.InterpolatedUnivariateSpline(e_tilde,
                                                                          interpolations,
                                                                          k=self._interpolation_degree)

            values = interpolator(energies)

        if self._nuFnu:

            return values / energies**2

        else:

            return values

    def to_dict(self, minimal=False):

        data = super(Function1D, self).to_dict(minimal)

        if not minimal:

            data['extra_setup'] = {'data_file': self._data_file}

        return data





