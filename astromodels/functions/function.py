import numpy as np

from astromodels.named_object import NamedObject

from astromodels.formula_parser import Formula

from astromodels.parameter import Parameter

import yaml

import parser

import collections


class FunctionDefinitionError(Exception):
    pass

def input_as_array(method):
    """
    Decorator which allows the decorated functions to be coded as having always a np.array as input. However, the
    decorated function will work with any input among a list, a single float or a np.array. For example::

      > def myfunc(x):
          idx = x < 1.0
          out = np.zeros_like(x)
          out[idx] = 0.0
          out[~idx] = 1.0
          return out

    This defines a function which returns 0 for all elements below 1 and 1 above that. As it is written it only works
    with a np.array as input. Indeed::

      > print( myfunc( np.array([0,2]) ) )
      [0,1]
      > print( myfunc(1.0) )
      IndexError: too many indices for array
      > print( myfunc([0,1,2,3]) )
      array([0, 0, 0, 1])

    Calling the function with a single float as input fails, while calling it with a list does not fail but returns the
    wrong result. If we decorate the function instead:

      > dec_func = input_always_array(myfunc)
      > print( dec_func(1.0) )
      1.0
      > print( dec_func([0,1,2,3]) )
      [0 1 1 1]
      > dec_func(np.array([0,2]))
      [0,1]

    :param method: method to decorate
    :return: same as the original method, but with the dimensions squeezed to the minimum. For example, an array with
    only one element will become a single number.
    """

    # Lookup in local scope is always much faster than in module scope. Hence, to decrease the performance hit of the
    # decorator, we store the functions we are going to use in the wrapper as local variables. Note that this part
    # of the code will be executed only once during the loading of the module using the decorator, while the wrapper
    # will be executed every time the method is called

    np_array = np.array
    np_squeeze = np.squeeze

    def wrapper(input_value, *args, **kwargs):

        # Transform the input in a numpy array, if needed.
        # If the input was a single float, this will become an array with shape
        # (1,), otherwise it will keep the shape of the input.

        if isinstance(input_value, np.ndarray):

            return method(input_value, *args, **kwargs)

        else:

            new_input = np_array(input_value, ndmin=1, copy=False)

            result = method(new_input, *args, **kwargs)

            # Now remove all dimensions of size 1. For example, an array of shape (1,) will become a single number.

            return np_squeeze(result)


    return wrapper

def clip(**myargs):
    """
    Force the output of the method to be always between min_value and max_value. For example, to constrain a function
    to be positive (i.e. all negative results will be substituted with 0):

      @clip(min_value=0)
      def myfun(in_value):
        return -1

    To constrain a function to be between 1e-20 and 1e20:
      @clip( min_value=1e-20, max_value=1e20 )
      def function(in):
        return -20 #will return 1e-20

    :param **myargs
    :keyword min_value: minimum allowed value (default: no boundary)
    :keyword max_value: maximum allowed value (default: no boundary)
    :return: same as the input method, but the array is clipped to min_value and max_value
    """

    min_value = None if 'min_value' not in myargs.keys() else myargs['min_value']
    max_value = None if 'max_value' not in myargs.keys() else myargs['max_value']

    def clip_generator(method):

        def wrapper(*args, **kwargs):

            result = method(*args, **kwargs)

            # This will clip the result array in-place

            result.clip(min_value,max_value,out=result)

            return result

        return wrapper

    return clip_generator

def output_always_finite(method):
    """
    This calls numpy.nan_to_num on the output of method, i.e., substitute nan with 0 and inf with a large number.

    :param method: method to decorate
    :return: same as the original method, but with nan transformed to 0 and inf to a large number.
    """

    np_nan_to_num = np.nan_to_num

    def wrapper(*args, **kwargs):

        result = method(*args, **kwargs)

        return np_nan_to_num( result )

    return wrapper


class DefinitionParser(object):

    def __init__(self, yaml_file):
        """

        :param yaml_file: definition of the function
        :return:
        """

        self.yaml_file = yaml_file

        # Parse the YAML file
        # (we don't catch any exception on purpose because if the file cannot be read there is nothing we can do)

        with open(self.yaml_file) as f:

            definitions = yaml.load(f)

        self.functions = []

        for key,value in definitions.iteritems():

            _thisParameters = {}

            for p, pval in value['init values'].iteritems():

                _thisParameters[p] = Parameter(p, pval)

            self.functions.append( SimpleFunction(value['description'],value['formula'], _thisParameters) )




class SimpleFunction(NamedObject):

    def __init__(self, description, formula, parameters):
        """

        :param description:
        :param formula:
        :param parameters:
        :return:
        """

        # Save description as documentation for this instance

        self.__doc__ = description

        # Parse the formula

        self._formula = Formula(formula)

        # Check that the parameters match between the definition received and the one needed by the formula

        set1 = set(parameters.keys())
        set2 = set(self._formula.parameters)

        if set1 != set2:

                # The parameters are different. Figure out who is missing and raise an exception accordingly

                if set1 > set2:

                    missing = set1 - set2

                    msg = "Parameters %s have init values but are not used in function" % ",".join(missing)

                else:

                    missing = set2 - set1

                    msg = "Parameters %s are used in function but do not have init values" % ",".join(missing)

                raise FunctionDefinitionError(msg)

