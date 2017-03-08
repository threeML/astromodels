import collections
import cPickle

from astromodels.utils.io import display
from astromodels.core.node_ctype import _Node
from astromodels.utils.valid_variable import is_valid_variable_name


class DuplicatedNode(Exception):
    pass


class ProtectedAttribute(RuntimeError):
    pass


class NonExistingAttribute(RuntimeWarning):
    pass


# This is necessary for pickle to be able to reconstruct a NewNode class (or derivate)
# during unpickling
class NewNodeUnpickler(object):

    def __call__(self, cls):

        instance = cls.__new__(cls)

        return instance


class Node(_Node):

    # This apparently dumb constructor is needed otherwise pickle will fail

    def __init__(self, name):

        if name == "units":

            import pdb;pdb.set_trace()

        assert is_valid_variable_name(name), "Illegal characters in name %s. You can only use letters and numbers, " \
                                             "and _" % name

        assert name != "name", "You cannot call a node 'name', it is reserved."

        _Node.__init__(self, name)

    # The next two methods are necessary for pickle to work

    def __reduce__(self):

        state = {}
        state['children'] = self._get_children()
        state['name'] = self.name
        state['__dict__'] = self.__dict__

        return NewNodeUnpickler(), (self.__class__,), state

    def __setstate__(self, state):

        # Set the children

        self._add_children(state['children'])

        # Set the name of this node

        self._change_name(state['name'])

        # Restore everything else

        for k in state['__dict__']:

            self.__dict__[k] = state['__dict__'][k]

    # This is necessary for copy.deepcopy to work
    def __deepcopy__(self, memodict={}):

        return cPickle.loads(cPickle.dumps(self))

    # This is used by dir() and by the autocompletion in Ipython
    def __dir__(self):

        # Get the names of the attributes of the class
        l = self.__class__.__dict__.keys()

        # Get all the children
        l.extend([child.name for child in self._get_children()])

        return l

    def to_dict(self, minimal=False):

        this_dict = collections.OrderedDict()

        for child in self._get_children():

            this_dict[child.name] = child.to_dict(minimal)

        return this_dict

    def _repr__base(self, rich_output):

        return "Node with name %s" % self.name

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


class OldNode(object):

    def __init__(self, name):

        self._children = collections.OrderedDict()
        self._parent = None

        assert is_valid_variable_name(name), "Illegal characters in name %s. You can only use letters and numbers, " \
                                             "and _" % name

        self._name = name

    @property
    def name(self):
        """
        Returns the name of the node

        :return: a string containing the name
        """
        return self._name

    @property
    def path(self):

        return ".".join(self._get_path())

    def _add_children(self, children):

        for child in children:

            self._add_child(child)

    def _add_child(self, new_child, name=None, add_attribute=True):

        new_child._set_parent(self)

        if name is None:

            name = new_child.name

        if name in self._children:

            raise DuplicatedNode("You cannot use the same name (%s) for different nodes" % name)

        self._children[name] = new_child

        if add_attribute:

            # Add also an attribute with the name of the new child, to allow access with a syntax like
            # node.child

            object.__setattr__(self, name, new_child)

    def _get_child(self, child_name):

        return self._children[child_name]

    def _remove_child(self, child_name):

        object.__delattr__(self, child_name)

        return self._children.pop(child_name)

    def _set_parent(self, parent):

        # We need to use directly the __setattr__ method because the normal self._children = ... will trigger
        # an exception, because of the __setattr__ method of the DualAccessClass which forbids changing
        # nodes. However, this might reassign the parent to a class

        object.__setattr__(self, "_parent", parent)

    def _get_parent(self):

        return self._parent

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

                this_child = this_child._get_child(key)

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
