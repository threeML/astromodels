__author__ = 'giacomov'

""" Implements a dictionary for which assignment can be lock, so that an operation like 'd[key] = value' will fail.
The main reason for this is the 'parameters' dictionary for the models. To assign a float 'value' to a parameter
the syntax to be used is parameters[parname].value = value and not parameters[parname] = value, which would overwrite
the parameter class with the float. """

import collections


class ParameterDictionary(collections.OrderedDict):

    __is_readonly = False

    def lock(self):

        self.__is_readonly = True

    def unlock(self):

        self.__is_readonly = False

    def __setitem__(self, *args, **kwargs):

        if not self.__is_readonly:

            super(ParameterDictionary, self).__setitem__(*args, **kwargs)

        else:

            raise ValueError("To assign values to the parameters use 'parameters[parname].value = value', and not "
                             "'parameters[parname] = value'")



