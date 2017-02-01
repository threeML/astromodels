import collections
import functools
import contextlib


_WITH_MEMOIZATION = True


@contextlib.contextmanager
def use_astromodels_memoization(switch):
    """
    Activate/deactivate memoization temporarily

    :param switch: True (memoization on) or False (memoization off)
    :return:
    """

    global _WITH_MEMOIZATION

    old_status = bool(_WITH_MEMOIZATION)

    _WITH_MEMOIZATION = bool(switch)

    yield

    _WITH_MEMOIZATION = old_status



def memoize(method):
    """
    A decorator for the 2d functions which memoize the results of 1 call (useful when the minimizer is taking partial
    derivatives and calls the 2d function several times with the same arguments)

    :param method: method to be memoized
    :return: the decorated method
    """

    cache = method.cache = collections.OrderedDict()

    # Put these two methods in the local space (faster)
    _get = cache.get
    _popitem = cache.popitem

    @functools.wraps(method)
    def memoizer(instance, x, *args, **kwargs):

        if not _WITH_MEMOIZATION:

            # Memoization is not active

            return method(instance, x, *args, **kwargs)

        # Create a tuple because a tuple is hashable

        unique_id = tuple(float(x.value) for x in instance.parameters.values()) + (x.size, x.min(), x.max())

        # Create a unique identifier for this combination of inputs

        key = hash(unique_id)

        # Let's do it this way so we only look into the dictionary once

        result = _get(key)

        if result is not None:

            return result

        else:

            result = method(instance, x, *args, **kwargs)

            cache[key] = result

            if len(cache) > 1000:
                # Remove the 100 elements that were put in first

                [_popitem(False) for i in range(100)]

            return result

    # Add the function as a "attribute" so we can access it
    memoizer.input_object = method

    return memoizer
