__author__ = 'giacomov'


class SpacesNotAllowedInName(ValueError):
    pass

class NamedObject(object):

    def __init__(self, name, allow_spaces=True):

        """

        :rtype : object
        """

        if not allow_spaces and name.find(" ") >= 0:

            raise SpacesNotAllowedInName("The name '%s' contains spaces, which are not allowed" % name)

        self.__name = name

    # Define the property "name" but make it read-only

    @property
    def name(self):
        """
        Get name of current object

        :return: the name of the object
        """
        return self.__name