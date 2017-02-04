import collections

import astropy.units as u
import numpy as np
import os
import pandas as pd
import re
import scipy.interpolate
import warnings
from pandas import HDFStore

from astromodels.core.parameter import Parameter
from astromodels.functions.function import Function1D, FunctionMeta
from astromodels.utils.configuration import get_user_data_path
from astromodels.core.my_yaml import my_yaml

# A very small number which will be substituted to zero during the construction
# of the templates
_TINY_ = 1e-50


__all__ = ["IncompleteGrid", "ValuesNotInGrid", "MissingDataFile", "TemplateModelFactory", "TemplateModel"]


class IncompleteGrid(RuntimeError):
    pass


class ValuesNotInGrid(ValueError):
    pass


class MissingDataFile(RuntimeError):
    pass


# This dictionary will keep track of the new classes already created in the current session
_classes_cache = {}


class TemplateModelFactory(object):

    def __init__(self, name, description, energies, names_of_parameters,
                 interpolation_degree=1, spline_smoothing_factor=0):

        # Store model name

        # Enforce that it does not contain spaces nor strange characters

        name = str(name)

        if re.match("[a-zA-Z_][a-zA-Z0-9_]*", name) is None:

            raise RuntimeError("The provided name '%s' is not a valid name. You cannot use spaces, "
                               "or special characters")

        self._name = name

        self._description = str(description)

        # Store energy grid

        if not isinstance(energies, u.Quantity):

            warnings.warn("Energy unit is not a Quantity instance, so units has not been provided. Using keV.")

            energies = energies * u.keV

        self._energies = np.array(energies.to(u.keV).value)

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

        self._spline_smoothing_factor = int(spline_smoothing_factor)

    def define_parameter_grid(self, parameter_name, grid):

        assert parameter_name in self._parameters_grids, "Parameter %s is not part of this model" % parameter_name

        grid_ = np.array(grid)

        assert grid_.shape[0] > 1, "A grid for a parameter must contain at least two elements"

        # Assert that elements are unique

        assert np.all(np.unique(grid_) == grid_), "Non-unique elements in grid for parameter %s" % parameter_name

        self._parameters_grids[parameter_name] = grid_

    def add_interpolation_data(self, differential_fluxes, **parameters_values_input):

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

        # Make sure we have all parameters and order the values in the same way as the dictionary
        parameters_values = np.zeros(len(self._parameters_grids)) * np.nan

        for key in parameters_values_input:

            assert key in self._parameters_grids, "Parameter %s is not known" % key

            idx = self._parameters_grids.keys().index(key)

            parameters_values[idx] = parameters_values_input[key]

        # If the user did not specify one of the parameters, then the parameters_values array will contain nan

        assert np.all(np.isfinite(parameters_values)), "You didn't specify all parameters' values."

        # Make sure we are dealing with pure numpy arrays (list and astropy.Quantity instances will be transformed)
        # First we transform the input into a u.Quantity (if it's not already)

        if not isinstance(differential_fluxes, u.Quantity):

            differential_fluxes = np.array(differential_fluxes) * 1 / (u.keV * u.s * u.cm ** 2)  # type: u.Quantity

        # Then we transform it in the right units and we cast it back to a pure np.array

        differential_fluxes = np.array(differential_fluxes.to(1 / (u.keV * u.s * u.cm ** 2)).value)

        # Now let's check for valid inputs

        assert self._energies.shape[0] == differential_fluxes.shape[0], "Differential fluxes and energies must have " \
                                                                        "the same number of elements"

        # Check that the provided value does not contains nan, inf nor zero (as the interpolation happens in the
        # log space)
        assert np.all(np.isfinite(differential_fluxes)), "You have invalid values in the differential flux (nan or inf)"
        assert np.all(differential_fluxes >= 0), "You have negative values in the differential flux (which is of " \
                                                 "course impossible)"

        if not np.all(differential_fluxes > 0):

            warnings.warn("You have zeros in the differential flux. Since the interpolation happens in the log space, "
                          "this cannot be accepted. We will substitute zeros with %g" % _TINY_)

            idx = (differential_fluxes == 0)  # type: np.ndarray
            differential_fluxes[idx] = _TINY_

        # Now set the corresponding values in the data frame

        # Now set the values in the data frame

        try:

            self._data_frame.loc[tuple(parameters_values)] = pd.to_numeric(differential_fluxes)

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
                                                             'spline_smoothing_factor': self._spline_smoothing_factor
                                                             }

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


def template_model_class_factory(model_name):
    """
    This factory creates a new class (not a new instance) corresponding to the given template. The new class can
    be used to instance a new model corresponding to the given template

    :param model_name: the name of the template
    :return: a new class called TemplateModel_[model_name]
    """

    # Create a new name for the new class
    new_name = 'TemplateModel_%s' % model_name

    # See if we already created this class
    if new_name in _classes_cache:

        return _classes_cache[new_name]

    else:

        # Create a new class for this template

        # Create an empty dictionary which will contain all methods/attributes of the new class
        new_class_dict = {}

        # Fill it with the methods/attributes of the TemplateModel class
        for item in TemplateModel.__dict__:

            if item == "__new__":

                # We don't want to copy the __new__ to avoid infinite recursion when creating the object

                continue

            new_class_dict[item] = TemplateModel.__dict__[item]

        # Now overwrite the things that specialize the class to this particular template

        # Read in template file
        # Get the data directory

        data_dir_path = get_user_data_path()

        # Sanitize the data file

        filename_sanitized = os.path.abspath(os.path.join(data_dir_path, '%s.h5' % model_name))

        if not os.path.exists(filename_sanitized):
            raise MissingDataFile("The data file %s does not exists. Did you use the "
                                  "TemplateFactory?" % (filename_sanitized))

        # Open the template definition and read from it

        new_class_dict['_data_file'] = filename_sanitized

        with HDFStore(filename_sanitized) as store:

            new_class_dict['_data_frame'] = store['data_frame']

            new_class_dict['_parameters_grids'] = collections.OrderedDict()

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

                    new_class_dict['_parameters_grids'][this_parameter_name] = store[key]

                    processed_parameters += 1

            new_class_dict['_energies'] = store['energies']

            # Now get the metadata

            metadata = store.get_storer('data_frame').attrs.metadata

            description = metadata['description']

            # name = metadata['name']

            new_class_dict['_interpolation_degree'] = metadata['interpolation_degree']

            new_class_dict['_spline_smoothing_factor'] = metadata['spline_smoothing_factor']

        # Make the dictionary of parameters

        function_definition = collections.OrderedDict()

        function_definition['description'] = description

        function_definition['latex'] = 'n.a.'

        # Now build the parameters according to the content of the parameter grid

        parameters_definition = collections.OrderedDict()

        parameters_definition['K'] = collections.OrderedDict((('desc', 'Normalization'), ('initial value', 1.0)))
        parameters_definition['scale'] = collections.OrderedDict((('desc', 'Normalization'),
                                                                ('initial value', 1.0),
                                                                ('min', 1e-5)))

        for parameter_name in new_class_dict['_parameters_grids'].keys():
            grid = new_class_dict['_parameters_grids'][parameter_name]

            parameters_definition[parameter_name] = collections.OrderedDict((('desc', 'none'),
                                                                           ('initial value', float(grid.median())),
                                                                           ('min', float(grid.min())),
                                                                           ('max', float(grid.max()))))

        function_definition['parameters'] = parameters_definition

        # Update the doc string so that the FunctionMeta class will populate correctly the class

        new_class_dict['__doc__'] = my_yaml.dump(function_definition)

        # Now susbsitute the evaluate function with a version with all the required parameters

        # Get the parameters' names (except for K and scale)
        par_names_no_K_no_scale = parameters_definition.keys()[2:]

        function_code = 'def new_evaluate(self, x, %s): ' \
                        'return K * self._interpolate(x, scale, [%s])' % (",".join(parameters_definition.keys()),
                                                                          ",".join(par_names_no_K_no_scale))

        exec (function_code)

        new_class_dict['_evaluate'] = new_evaluate

        new_class_dict['evaluate'] = new_evaluate

        new_class = FunctionMeta(new_name, (Function1D,), new_class_dict)

        # Now prepare the interpolators
        new_class._prepare_interpolators()

        return new_class


class TemplateModelUnplickler(object):

    def __call__(self, model_name):

        new_instance = TemplateModel(model_name)

        return new_instance


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
                       at E = scale * x.This is useful for example when the template describe energies in the rest
                       frame, at which point the scale describe the transformation between rest frame energy and
                       observer frame energy. Fix this to 1 to neutralize its effect.

                initial value : 1.0
                min : 1e-5

        """

    # Use __new__ instead of init so we can use the factory
    def __new__(cls, model_name):

        new_class = template_model_class_factory(model_name)

        new_instance = new_class()

        # Add attributes from the type to the instance (for easier accessibility)
        new_instance._data_file = new_class._data_file

        new_instance._data_frame = new_class._data_frame

        new_instance._parameters_grids = new_class._parameters_grids

        new_instance._energies = new_class._energies

        new_instance._interpolation_degree = new_class._interpolation_degree

        new_instance._spline_smoothing_factor = new_class._spline_smoothing_factor

        new_instance._interpolators = new_class._interpolators

        new_instance._model_name = model_name

        return new_instance

    # This is a class method because we only need to make this once for a given template class
    # (no need to repeat this for every instance)
    @classmethod
    def _prepare_interpolators(cls):

        # Figure out the shape of the data matrices
        data_shape = map(lambda x: x.shape[0], cls._parameters_grids.values())

        cls._interpolators = []

        for energy in cls._energies:

            # Make interpolator for this energy
            # NOTE: we interpolate on the logarithm

            this_data = np.array(np.log10(cls._data_frame[energy].values).reshape(*data_shape), dtype=float)

            if len(cls._parameters_grids.values()) == 2:

                x, y = cls._parameters_grids.values()

                # Make sure that the requested polynomial degree is less than the number of data sets in
                # both directions

                msg = "You cannot use an interpolation degree of %s if you don't provide at least %s points " \
                      "in the %s direction. Increase the number of templates or decrease the interpolation " \
                      "degree."

                if len(x) <= cls._interpolation_degree:

                    raise RuntimeError(msg % (cls._interpolation_degree, cls._interpolation_degree + 1, 'x'))

                if len(y) <= cls._interpolation_degree:

                    raise RuntimeError(msg % (cls._interpolation_degree, cls._interpolation_degree + 1, 'y'))

                this_interpolator = RectBivariateSplineWrapper(x, y, this_data,
                                                               kx=cls._interpolation_degree,
                                                               ky=cls._interpolation_degree,
                                                               s=cls._spline_smoothing_factor)

            else:

                # In more than 2d we can only use linear interpolation

                this_interpolator = scipy.interpolate.RegularGridInterpolator(cls._parameters_grids.values(),
                                                                              this_data)

            cls._interpolators.append(this_interpolator)

    def _set_units(self, x_unit, y_unit):

        self.K.unit = y_unit

        self.scale.unit = 1 / x_unit

        # The other parameters must be scale free (the results of the interpolation must be unit free)

    # This function will be substituted during construction by another version with
    # all the parameters of this template

    def evaluate(self, x, K, scale):

        # This is overridden in the constructor

        raise NotImplementedError("Should not get here!")

    def _interpolate(self, energies, scale, parameters_values):

        if isinstance(energies, u.Quantity):

            # Templates are always saved with energy in keV. We need to transform it to
            # a dimensionless quantity (actually we take the .value property) because otherwise
            # the logarithm below will fail.

            energies = np.array(energies.to('keV').value, ndmin=1, copy=False, dtype=float)

            # Same for the scale

            scale = scale.to(1 / u.keV).value

        log_energies = np.log10(energies)

        e_tilde = self._energies * scale

        # Gather all interpolations for these parameters' values at all defined energies
        # (these are the logarithm of the values)

        log_interpolations = np.array(map(lambda i:self._interpolators[i](parameters_values),
                                          range(self._energies.shape[0])))

        # Now interpolate the interpolations to get the flux at the requested energies

        # NOTE: the variable "interpolations" contains already the log10 of the values,

        interpolator = scipy.interpolate.InterpolatedUnivariateSpline(np.log10(e_tilde),
                                                                      log_interpolations,
                                                                      k=self._interpolation_degree,
                                                                      ext=0)

        values = np.power(10, interpolator(log_energies))

        # The division by scale results from the differential:
        # E = e * scale
        # de = dE / scale
        # dN / dE = dN / de * de / dE = dN / de * (1 / scale)

        # NOTE: the units are added back through the multiplication by K in the evaluate method

        return values / scale

    @property
    def data_file(self):

        return self._data_file

    def to_dict(self, minimal=False):

        data = super(Function1D, self).to_dict(minimal)

        return data

    # Now implement custom methods for pickling
    def __reduce__(self):

        model_name = self._name.replace("TemplateModel_","")

        return TemplateModelUnplickler(), (model_name,), self.parameters

    def __setstate__(self, state):

        # Restore parameters
        for parameter_name in state:

            self.parameters[parameter_name] = state[parameter_name]
