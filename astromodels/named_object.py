__author__ = 'giacomov'

class NamedObject(object):

    def __init__(self, name):

        self._name = name

    # Define the property "name" but make it read-only

    @property
    def name(self):
        """
        Get name of current object

        :return: the name of the object
        """
        return self.__name