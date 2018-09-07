from builtins import range
import collections
import functools
import contextlib
import astropy.units as u

_WITH_MEMOIZATION = False
_CACHE_SIZE = 20


@contextlib.contextmanager
def use_astromodels_memoization(switch, cache_size=_CACHE_SIZE):
    """
    Activate/deactivate memoization temporarily

    :param switch: True (memoization on) or False (memoization off)
    :param cache_size: number of previous evaluation of functions to keep in memory. Default: 100
    :return:
    """

    global _WITH_MEMOIZATION
    global _CACHE_SIZE

    old_status = bool(_WITH_MEMOIZATION)
    old_cache_size = int(_CACHE_SIZE)

    _WITH_MEMOIZATION = bool(switch)
    _CACHE_SIZE = int(cache_size)

    yield

    _WITH_MEMOIZATION = old_status
    _CACHE_SIZE = old_cache_size



def memoize(method):
    """
    A decorator for functions of sources which memoize the results of the last _CACHE_SIZE calls,

    :param method: method to be memoized
    :return: the decorated method
    """

    cache = method.cache = collections.OrderedDict()

    # Put these two methods in the local space (faster)
    _get = cache.get
    _popitem = cache.popitem

    @functools.wraps(method)
    def memoizer(instance, x, *args, **kwargs):

        if not _WITH_MEMOIZATION or isinstance(x, u.Quantity):

            # Memoization is not active or using units, do not use memoization

            return method(instance, x, *args, **kwargs)

        # Create a tuple because a tuple is hashable

        unique_id = tuple(float(yy.value) for yy in list(instance.parameters.values())) + (x.size, x.min(), x.max())

        # Create a unique identifier for this combination of inputs

        key = hash(unique_id)

        # Let's do it this way so we only look into the dictionary once

        result = _get(key)

        if result is not None:

            return result

        else:

            result = method(instance, x, *args, **kwargs)

            cache[key] = result

            if len(cache) > _CACHE_SIZE:
                # Remove half of the element (but at least 1, even if _CACHE_SIZE=1, which would be pretty idiotic ;-) )
                [_popitem(False) for i in range(max(_CACHE_SIZE // 2, 1))]

            return result

    # Add the function as a "attribute" so we can access it
    memoizer.input_object = method

    return memoizer
