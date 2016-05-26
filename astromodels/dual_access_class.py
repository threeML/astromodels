__author__ = 'giacomov'

from astropy.units import Quantity
import warnings


class ProtectedAttribute(RuntimeError):
    pass


class NonExistingAttribute(RuntimeWarning):
    pass


class DualAccessClass(object):
    """
    Suppose there is a class A having a dictionary b as a member. Inheriting from this class will allow the user to
    access the elements of b either as A.b[key] or as A.key. Accessing the element as A.key is read-only,
    i.e., A.key=something will fail.
    """

    def __init__(self, label, dictionary):

        self._lookup_dictionary = dictionary
        self._dictionary_label = label

        for key, value in self._lookup_dictionary.iteritems():

            super(DualAccessClass, self).__setattr__(key, value)

    def __setattr__(self, key, value):

        if hasattr(self, '_lookup_dictionary') and key in self._lookup_dictionary:

            item = self._lookup_dictionary[key]

            if hasattr(item, 'value') and not isinstance(item, Quantity):

                item.value = value

            else:

                raise ProtectedAttribute("You cannot assign to a %s" % self._dictionary_label)

        else:

            # Attributes which start with "_" are created by children classes

            if not key[0] == '_' and not hasattr(self, key):

                warnings.warn("Attribute %s does not exist. Check for typos." % key, NonExistingAttribute)

            super(DualAccessClass, self).__setattr__(key, value)

    def _add_attribute(self, name, value):

        self._lookup_dictionary[name] = value

        super(DualAccessClass, self).__setattr__(name, value)

    def _del_attribute(self, name):

        super(DualAccessClass, self).__delattr__(name)

        return self._lookup_dictionary.pop(name)
