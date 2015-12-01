import numpy as np

def input_always_array(method):
    """
    Decorator which allows to always handle the input as a numpy.array, even if the function is called with a single
    float as argument. A function decorated with this can be called in the same way with the following inputs: one
    single number, a python iterable like a list, or a numpy.array. For example::

      > myfunc = lambda x: x * 2.0
      > dec_func = input_always_array(myfunc)
      > print( dec_func(1) )
      2.0
      > print( dec_func([1,2]) )
      array(

    :param method: method to decorate
    :return: same as the original method, but with the dimensions squeezed to the minimum. For example, an array with
    only one element will become a single number.
    """

    def wrapper(input_value, *args, **kwargs):

        # Transform the input in a numpy array. If the input was a single float, this will become an array with shape
        # (1,), otherwise it will keep the shape of the input. Note that the performance hit is very minimal, of the
        # order of a few microseconds

        new_input = np.array(input_value, ndmin=1)

        result = method(new_input, *args, **kwargs)

        # Now remove all dimensions of size 1. For example, an array of shape (1,) will become a single number.

        return np.squeeze(result)

    return wrapper

def clip(**myargs):
    """
    Force the output of the method to be always between min_value and max_value. For example, to constrain a function
    to be positive (i.e. all negative results will be substituted with 0):

      @clip(min_value=0)
      def function(in):
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

    def wrapper(*args, **kwargs):

        result = method(*args, **kwargs)

        return numpy.nan_to_num( result )

    return wrapper
