from builtins import object
import numpy as np


class ParameterTransformation(object):
    def __init__(self, is_positive=False):

        self._is_positive = is_positive

    @property
    def is_positive(self):
        return self._is_positive
    
    def forward(self, external_value):

        raise NotImplementedError("You have to implement this")

    def backward(self, internal_value):

        raise NotImplementedError("You have to implement this")


class LogarithmicTransformation(ParameterTransformation):

    def __init__(self):

        super(LogarithmicTransformation, self).__init__(is_positive=True)
    
    def forward(self, external_value):

        #  Throw an error if taking the logarithm of a negative number (or nan)

        with np.errstate(invalid='raise'):

            res = np.log10(external_value)

        return res

    def backward(self, internal_value):

        return 10**internal_value


_known_transformations = {'log10': LogarithmicTransformation}


def get_transformation(transformation_name):
    """
    Returns an instance of a transformation by name

    :param transformation_name:
    :return: instance of transformation with provided name
    """


    if not transformation_name in _known_transformations:

        raise ValueError("Transformation %s is not known" % transformation_name)

    else:

        return _known_transformations[transformation_name]()
