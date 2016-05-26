__author__ = 'giacomov'

import collections

from astromodels.utils.io import display
from astromodels.dual_access_class import DualAccessClass


class DuplicatedNode(Exception):
    pass


def name_is_legal(name):
    """
    Check that the name is legal, i.e., it does not contain special characters nor spaces

    :param name: the name to check
    :return: True or False
    """

    # Define the set of valid characters
    valid_char_sequence = 'qwertyuiopasdfghjklzxcvbnmQWERTYUIOPASDFGHJKLZXCVBNM1234567890-+_'

    valid = set(valid_char_sequence).issuperset(name)

    if valid:

        # Check that the name does not start with -
        if name[0] == '-':

            return False

    return valid


class Node(DualAccessClass):

    def __init__(self, name):

        self.__children = collections.OrderedDict()
        self.__parent = None

        assert name_is_legal(name),"Illegal characters in name %s. You can only use letters and numbers, " \
                                   "and + and - (but the name cannot start with -)" % name

        self._name = name

        super(Node, self).__init__('node', self.__children)

    @property
    def name(self):
        """
        Returns the name of the node

        :return: a string containing the name
        """
        return self._name

    @property
    def _children(self):
        return self.__children

    @property
    def path(self):

        return ".".join(self._get_path())

    def _reset_node(self):

        self.__children = collections.OrderedDict()
        self.__parent = None

    def _add_children(self, children):

        for child in children:

            self._add_child(child)

    def _add_child(self, new_child, name=None):

        new_child._set_parent(self)

        if name is None:

            name = new_child.name

        if name in self.__children:

            raise DuplicatedNode("You cannot use the same name (%s) for different nodes" % name)

        self.__children[name] = new_child

        # Add also an attribute with the name of the new child, to allow access with a syntax like
        # node.child

        self._add_attribute(name, new_child)

    def _get_child(self, child_name):

        return self.__children[child_name]

    def _remove_child(self, child_name):

        return self._del_attribute(child_name)

    def _set_parent(self, parent):

        self.__parent = parent

    def _get_parent(self):

        return self.__parent

    def _get_child_from_path(self, path):
        """
        Return a children below this level, starting from a path of the kind "this_level.something.something.name"

        :param path: the key
        :return: the child
        """

        keys = path.split(".")

        this_child = self

        for key in keys:

            try:

                this_child = this_child.get_child(key)

            except KeyError:

                raise KeyError("Child %s not found" % path)

        return this_child

    def _get_path(self):

        parent_names = []

        current_node = self

        while True:

            this_parent = current_node._get_parent()

            if this_parent is None:

                break

            else:

                parent_names.append(current_node.name)

                current_node = this_parent

        return parent_names[::-1]

    def to_dict(self, minimal=False):

        this_dict = collections.OrderedDict()

        for key, val in self.__children.iteritems():

            this_dict[key] = val.to_dict(minimal)

        return this_dict

    def _repr__base(self, rich_output):

        raise NotImplementedError("You should implement the __repr__base method for each class")

    def __repr__(self):
        """
        Textual representation for console

        :return: representation
        """

        return self._repr__base(rich_output=False)

    def _repr_html_(self):
        """
        HTML representation for the IPython notebook

        :return: HTML representation
        """

        return self._repr__base(rich_output=True)

    def display(self):
        """
        Display information about the point source.

        :return: (none)
        """

        # This will automatically choose the best representation among repr and repr_html

        display(self)

    def _find_instances(self, cls):
        """
        Find all the instances of cls below this node.

        :return: a dictionary of instances of cls
        """

        instances = collections.OrderedDict()

        for child_name, child in self._children.iteritems():

            if isinstance(child, cls):

                key_name = ".".join(child._get_path())

                instances[key_name] = child

                # Now check if the instance has children,
                # and if it does go deeper in the tree

                # NOTE: an empty dictionary evaluate as False

                if child._children:

                    instances.update(child._find_instances(cls))

            else:

                instances.update(child._find_instances(cls))

        return instances
