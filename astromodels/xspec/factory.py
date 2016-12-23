import collections
import sys

import astropy.units as u
import os
import re
import warnings

from astromodels.core.my_yaml import my_yaml
from astromodels.functions.function import get_function_class
from astromodels.utils.configuration import get_user_data_path


class XSpecNotAvailable(ImportWarning):
    pass

try:

    from astromodels.xspec import _xspec

except ImportError:

    raise XSpecNotAvailable("You need to have XSPEC installed and configured in order to use its models.")


# This list defines all python protected names (names which variables or attributes should not have)
illegal_variable_names_ = 'and, assert, break, class, continue, def, del, elif, else, except, exec, finally, for,' \
                         'if, import, in, is, lambda, not, or, pass, print, raise, return, try, while, data, ' \
                         'float, int, numeric, array, close, input, open, range, type, write, zeros, from, global'
# Make it a list
illegal_variable_names = illegal_variable_names_.replace(" ","").split(",")


def find_model_dat():
    """
    Find the file containing the definition of all the models in Xspec
    (model.dat) and return its path
    """

    # model.dat is in $HEADAS/../spectral

    headas_env = os.environ.get("HEADAS")

    assert headas_env is not None, ("You need to setup the HEADAS variable before importing this module."
                                    " See Heasoft documentation.")

    # Expand all variables and other things like ~
    headas_env = os.path.expandvars(os.path.expanduser(headas_env))

    # Lazy check that it exists

    assert os.path.exists(headas_env), "The HEADAS env. variable point to a non-existent directory: %s" % (headas_env)

    # Get one directory above HEADAS (i.e., $HEADAS/..)

    inferred_path = os.path.dirname(headas_env)

    # Now model.dat should be in $HEADAS/../spectral/manager

    final_path = os.path.join(inferred_path, 'spectral', 'manager', 'model.dat')

    # Check that model.dat exists

    assert os.path.exists(final_path), "Cannot find Xspec model definition file %s" % (final_path)

    return os.path.abspath(final_path)


def get_models(model_dat_path):
    """
    Parse the model.dat file from Xspec and returns a dictionary containing the definition of all the models

    :param model_dat_path: the path to the model.dat file
    :return: dictionary containing the definition of all XSpec models
    """

    # Check first if we already have a model data file in the data directory


    with open(model_dat_path) as f:

        # model.dat is a text file, no size issues here (will fit in memory)

        model_dat = f.read()

    # Replace \r or \r\n with \n (the former two might be in the file if
    # it originates from Windows)
    if "\n" in model_dat:

        model_dat = model_dat.replace("\r", "")

    else:

        model_dat = model_dat.replace("\r", "\n")

    # Loop through the file and build the model definition
    # dictionary

    lines = model_dat.split("\n")

    model_definitions = collections.OrderedDict()

    for line in lines:

        match = re.match('''(.+(add|mul|con|acn).+)''', line)

        if match is not None:

            # This is a model definition

            tokens = line.split()

            if len(tokens) == 7:

                (model_name, n_parameters,
                 min_energy, max_energy,
                 library_function,
                 model_type,
                 flag) = line.split()

            else:

                (model_name, n_parameters,
                 min_energy, max_energy,
                 library_function,
                 model_type,
                 flag,
                 flag_2) = line.split()

            this_model = collections.OrderedDict()

            this_model['description'] = 'The %s model from XSpec (https://heasarc.gsfc.nasa.gov/xanadu/' \
                                                   'xspec/manual/XspecModels.html)' % model_name

            this_model['parameters'] = collections.OrderedDict()

            model_definitions[(model_name, library_function, model_type)] = this_model

        else:

            # This is a parameter definition

            if len(line.split()) == 0:
                # Empty line
                continue

            # Parameters are free by default, unless they have delta < 0
            # they are switch parameters or scale parameters

            free = True

            if line[0] == '$':

                # Probably a switch parameter

                free = False

                tokens = line.split()

                if len(tokens) == 2:

                    par_name = tokens[0][1:]

                    default_value = tokens[1]

                    par_unit = ""

                    hard_minimum, soft_minimum, soft_maximum, hard_maximum = (0, 0, 1e9, 1e9)

                elif len(tokens) == 3:

                    par_name = tokens[0][1:]

                    default_value = tokens[2]

                    par_unit = ""

                    hard_minimum, soft_minimum, soft_maximum, hard_maximum = (0, 0, 1e9, 1e9)

                else:

                    match = re.match('(\S+)\s+(\".+\"|[a-zA-Z]+)?(.+)*', line[1:])

                    par_name, par_unit, par_spec = match.groups()

                    tokens = par_spec.split()

                    if len(tokens) == 1:

                        default_value = tokens[0]
                        par_unit = ""
                        hard_minimum, soft_minimum, soft_maximum, hard_maximum = (0, 0, 1e9, 1e9)

                    else:

                        par_unit = ""

                        (default_value,
                         hard_minimum, soft_minimum,
                         soft_maximum, hard_maximum,
                         delta) = par_spec.split()


            else:

                # A normal parameter

                match = re.match('(\S+)\s+(\".+\"|\S+)(.+)', line)

                if match is None:

                    raise RuntimeError("Cannot parse parameter %s" % line)

                par_name, par_unit, par_spec = match.groups()

                if par_name[0] == '*':

                    # Scale parameter (always frozen)

                    par_name = par_name[1:]

                    free = False

                    tokens = par_spec.split()

                    if len(tokens) == 1:

                        default_value = tokens[0]
                        (hard_minimum, soft_minimum,
                         soft_maximum, hard_maximum,
                         delta) = (None, None, None, None, 0.1)

                    else:

                        (default_value,
                         hard_minimum, soft_minimum,
                         soft_maximum, hard_maximum,
                         delta) = par_spec.split()

                        delta = abs(float(delta))

                else:

                    (default_value,
                     hard_minimum, soft_minimum,
                     soft_maximum, hard_maximum,
                     delta) = par_spec.split()

                    delta = float(delta)

                    if delta <= 0:
                        free = False

                        delta = abs(delta)

            # Now fix the parameter name removing illegal characters
            # (for example 'log(A)' is not a legal name)

            par_name = re.sub('[\(,\)]', '_', par_name)
            par_name = re.sub('<', '_less_', par_name)
            par_name = re.sub('>', '_more_', par_name)
            par_name = re.sub('/', '_div_', par_name)
            par_name = re.sub('\-', '_minus_', par_name)
            par_name = re.sub('\+', '_plus_', par_name)
            par_name = re.sub('\.', '_dot_', par_name)
            par_name = re.sub('@', '_at_', par_name)

            # Parameter names must be lower case
            par_name = par_name.lower()

            # Some parameters are enclosed between ", like "z"
            par_name = par_name.replace('"', '')

            # "z" is a protected name in astromodels.
            if par_name == "z":

                par_name = "redshift"

            # Check that the parameter name is not an illegal Python name
            if par_name in illegal_variable_names:

                par_name = "xs_%s" % par_name

            # Sometimes the unit is " " which is not recognized by astropy
            if par_unit:

                par_unit = par_unit.replace("\"",'')

            # There are some errors in model.dat , like KeV instead of keV, and ergcm/s instead of erg cm /s
            # Let's correct them
            if par_unit == "KeV":

                par_unit = "keV"

            elif par_unit == "ergcm/s":

                par_unit = "erg cm / s"

            elif par_unit == "days":

                par_unit = "day"

            elif par_unit == "z" and par_name.lower() == "redshift":

                par_unit = ""

            # There are funny units in model.dat, like "Rs" (which means Schwarzschild radius) or other things
            # so let's try to convert the par_unit into an astropy.Unit instance. If that fails, use a unitless unit
            try:

                _ = u.Unit(par_unit)

            except ValueError:

                # Unit not recognized

                #warnings.warn("Unit %s is not recognized by astropy." % par_unit)

                par_unit = ''

            # Make sure that this is a valid python identifier
            # by matching it with the relative regexp
            if re.match('([a-zA-Z_][a-zA-Z0-9_]*)$', par_name) is None:

                raise ValueError("Illegal identifier name %s" % (par_name))

            this_model['parameters'][par_name] = {'initial value': float(default_value),
                                                  'desc': '(see https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/'
                                                          'XspecModels.html)',
                                                  'min': hard_minimum,
                                                  'max': hard_maximum,
                                                  'delta': float(delta),
                                                  'unit': par_unit,
                                                  'free': free}

    return model_definitions

class_definition_code = '''

from astromodels.functions.function import FunctionMeta, Function1D
import numpy as np
import astropy.units as u
from astromodels.xspec import _xspec

class XS_$MODEL_NAME$(Function1D):

    """
$DOCSTRING$
    """

    __metaclass__ = FunctionMeta

    def _setup(self):

        # Link to the Xspec function
        self._model = _xspec.$XSPEC_FUNCTION$
        self._model_type = '$MODEL_TYPE$'

        if self._model_type == 'add':

            self._scale = 1e6

            self._fixed_units = (u.keV, 1 / (u.keV * u.cm**2 * u.s))

        else:

            # For multiplicative models, there is no differentiation to be made

            self._scale = 10

            self._fixed_units = (u.keV, u.dimensionless_unscaled)

    def evaluate(self, x, $PARAMETERS_NAMES$):

        quantity = False

        if isinstance(x, u.Quantity):

            x = np.array(x.to('keV').value, ndmin=1, copy=False, dtype=float)

            quantity = True

            parameters_tuple = tuple(map(lambda x:x.value, ($PARAMETERS_NAMES$,)))

        else:

            # Create a tuple of the current values of the parameters

            parameters_tuple = ($PARAMETERS_NAMES$,)

        # Finite difference differentiation of the Xspec
        # function for additive model. Indeed, xspec function return the integral
        # of the function on energy ranges, while we need the
        # differential flux. For multiplicative models instead Xspec returns just
        # the value for the multiplicative factor at the average between emin and
        # emax

        # Adapt the epsilon to the value to reduce the error

        epsilon = x / self._scale

        # In the normal case, when x is an array of more than one
        # element, the first call will succeed. If however there
        # is only one element in the array xspec will complain,
        # so handle that as a special case

        try:

            val = self._model(parameters_tuple, x - epsilon, x + epsilon)

        except TypeError:

            assert x.shape[0]==1, "This is a bug, xspec call failed and x is not only one element"
            
            val = self._model(parameters_tuple, ((x - epsilon)[0], (x + epsilon)[0]))[0]

        # val is now F(x-epsilon,x+epsilon) ~ f(x) * ( 2 * epsilon )

        if self._model_type == 'add':

            # In a additive model the function returns the integral over the bins

            final_value = val / (2 * epsilon)

        else:

            # In a multiplicative model the function returns the average factor over
            # the bins

            final_value = val

        if quantity:

            if self._model_type == 'add':

                return final_value * u.Unit('1 / (keV cm^2 s)')

            else:

                return final_value * u.dimensionless_unscaled

        else:

            return final_value

    def has_fixed_units(self):

        return True

    def _set_units(self, x_unit, y_unit):

        # Make sure this is an energy
        assert str(x_unit.physical_type) == 'energy', "Trying to set x-unit with (%s), which is not energy" % (x_unit)

        # Make sure the y_unit is the correct one
        try:

            y_unit.in_units(self._fixed_units[1])

        except:

            raise RuntimeError("Xspec model %s cannot have units of %s" % (self.name, y_unit))

    def _integral(self, low_bounds, hi_bounds, $PARAMETERS_NAMES$):

        # Create a tuple of the current values of the parameters
        parameters_tuple = ($PARAMETERS_NAMES$,)

        return self._model(parameters_tuple, low_bounds, hi_bounds)

'''


def xspec_model_factory(model_name, xspec_function, model_type, definition):

    class_name = 'XS_%s' % model_name

    # Get the path to the user data directory
    user_data_path = get_user_data_path()

    # Check if the code for this function already exists

    code_file_name = os.path.join(user_data_path, '%s.py' % class_name)

    if os.path.exists(code_file_name):

        # Code already exists
        pass

    else:

        print("Generating code for Xspec model %s..." % model_name)

        # If this is an additive model (model_type == 'add') we need to add
        # one more parameter (normalization)

        if model_type == 'add':

            definition['parameters']['norm'] = {'initial value': 1.0,
                                                      'desc': '(see https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/'
                                                              'XspecModels.html)',
                                                      'min': 0,
                                                      'max': None,
                                                      'delta': 0.1,
                                                      'unit': 'keV / (cm2 s)',
                                                      'free': True}

        assert model_type != 'con', "Convolution models are not yet supported"

        # Get a list of the parameter names
        parameters_names = ", ".join(definition['parameters'].keys())

        # Create the docstring
        docstring = my_yaml.dump(definition)

        # Create the class by substituting in the class_definition_code the
        # relevant things for this model

        code = class_definition_code.replace('$MODEL_NAME$', model_name)
        code = code.replace('$DOCSTRING$', docstring)
        code = code.replace('$PARAMETERS_NAMES$', parameters_names)
        code = code.replace('$XSPEC_FUNCTION$', xspec_function)
        code = code.replace('$MODEL_TYPE$', model_type)

        # Write to the file

        with open(code_file_name, 'w+') as f:

            f.write("# This code has been automatically generated. Do not edit.\n")
            f.write("\n\n%s\n" % code)

    # Add the path to sys.path if it doesn't
    if user_data_path not in sys.path:

        sys.path.append(user_data_path)

    # Import the class in the current namespace (locals)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        exec('from %s import %s' % (class_name, class_name))

    # Return the class we just created

    return class_name, locals()[class_name]


def setup_xspec_models():

    classes = []

    sys.stdout.write("Loading xspec models...")

    all_models = get_models(find_model_dat())

    for (model_name, xspec_function, model_type) in all_models:

        if model_type == 'con':

            # convolution models are not supported
            continue

        if not hasattr(_xspec, xspec_function):

            # Some function do not exist in the wrapper. Let's ignore them

            continue

        this_model = all_models[(model_name, xspec_function, model_type)]

        # When the class is created it is registered among the known functions in the function module
        # (it happens in the metaclass), so we don't need to do anything special here after the
        # class type is created

        this_class_name, this_class = xspec_model_factory(model_name, xspec_function, model_type, this_model)

        classes.append(this_class_name)

    sys.stdout.write("done\n")

    return classes

# This will either work or issue a warning if XSpec is not available

new_functions = setup_xspec_models()

# Now import the new classes in the local namespace (if any)
# This is needed to make the classes pickeable

__all__ = []

for function_name in new_functions:

    __all__.append(function_name)

    locals()[function_name] = get_function_class(function_name)
