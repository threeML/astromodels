__author__ = 'giacomov'


class ProtectedAttribute(RuntimeError):
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

        for key,value in self._lookup_dictionary.iteritems():

            super(DualAccessClass, self).__setattr__(key, value)

    def __setattr__(self, key, value):

        if hasattr(self, '_lookup_dictionary') and key in self._lookup_dictionary:

            raise ProtectedAttribute("You cannot assign to a %s" % self._dictionary_label)

        else:

            super(DualAccessClass, self).__setattr__(key, value)