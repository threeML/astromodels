__author__ = 'giacomov'

import collections

from astromodels.utils.io import display
from astromodels.dual_access_class import DualAccessClass


class DuplicatedNode(Exception):
    pass


class Node(DualAccessClass):

    def __init__(self):

        self._children = collections.OrderedDict()
        self._parent = None
        
        super(Node, self).__init__('node', self._children)

    @property
    def children(self):
        return self._children

    @property
    def path(self):

        return ".".join(self.get_path())

    def reset(self):

        self._children = collections.OrderedDict()
        self._parent = None

    def add_children(self, children):

        for child in children:

            self.add_child(child)

    def add_child(self, new_child, name=None):

        new_child.set_parent(self)

        if name is None:

            name = new_child.name

        if name in self._children:

            raise DuplicatedNode("You cannot use the same name (%s) for different nodes" % name)

        self._children[name] = new_child

        # Add also an attribute with the name of the new child, to allow access with a syntax like
        # node.child

        self.add_attribute(name, new_child)

    def get_child(self, child_name):

        return self._children[child_name]

    def remove_child(self, child_name):

        return self.del_attribute(child_name)

    def get_children(self):

        return self._children

    def set_parent(self, parent):

        self._parent = parent

    def get_parent(self):

        return self._parent

    def get_child_from_path(self, path):
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

    def get_path(self):

        parent_names = []

        current_node = self

        while True:

            this_parent = current_node.get_parent()

            if this_parent is None:

                break

            else:

                parent_names.append(current_node.name)

                current_node = this_parent

        return parent_names[::-1]

    def to_dict(self, minimal=False):

        this_dict = collections.OrderedDict()

        for key, val in self._children.iteritems():

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

    def find_instances(self, cls):
        """
        Find all the parameters below this node. Note that given the flexible way in which parameters can be defined,
        linked and so on, they can be in many different places in the tree.

        :return: a dictionary of parameters
        """

        # NOTE: parameters are always at the end of a branch in the tree. We exploit this fact
        # so we can avoid using isinstance() or similar and we keep duck typing

        parameters = collections.OrderedDict()

        for child_name, child in self.get_children().iteritems():

            if isinstance(child, cls):

                key_name = ".".join(child.get_path())

                parameters[key_name] = child

                # Now check if the parameter has children,
                # and if it does go deeper in the tree

                # NOTE: an empty dictionary evaluate as False

                if child.get_children():

                    parameters.update(child.find_instances(cls))

            else:

                parameters.update(child.find_instances(cls))

        return parameters
