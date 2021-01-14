import collections
import itertools
from builtins import object
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Type

from future import standard_library

from astromodels.utils.io import display
from astromodels.utils.valid_variable import is_valid_variable_name

from .cpickle_compatibility_layer import cPickle

standard_library.install_aliases()


#from astromodels.core.node_ctype import _Node


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


@dataclass(repr=False)
class _Node:
    _name: str
    _parent: Optional[Type["_Node"]] = field(repr=False, default=None)
    _children: Dict[str, Type["_Node"]] = field(
        default_factory=dict, repr=False)
    _path: Optional[str] = field(repr=False, default="")

    # The next 3 methods are *really* necessary for anything to work

    def __reduce__(self):

        state = {}
        state['parent'] = self._get_parent()
        state['path'] = self._path
        state['children'] = self._get_children()
        state['name'] = self._name
        state['__dict__'] = self.__dict__

        return NewNodeUnpickler(), (self.__class__,), state

    def __setstate__(self, state):

        self._children = {}

        # Set the name of this node

        self._name = state['name']

        # Set the parent

        self._parent = state['parent']

        # set the path

        self._path = state['path']

        # Set the children

        self._add_children(state['children'])

        # Restore everything else

        for k in state['__dict__']:

            self.__dict__[k] = state['__dict__'][k]

    # This is necessary for copy.deepcopy to work
    def __deepcopy__(self, memodict={}):

        return cPickle.loads(cPickle.dumps(self))

    def _add_child(self, child: Type["_Node"]):

        try:
            if child._name not in self._children:
                self._children[child._name] = child
                child._set_parent(self)

            else:

                raise RuntimeError(
                    f"A child with name {child._name} already exists")

        except(AttributeError):

            raise TypeError()

    def _add_children(self, children: List[Type["_Node"]]):

        for child in children:
            self._add_child(child)

    def _remove_child(self, name: str):
        """
        return a child
        """
        del self._children[name]
        # return self._children.pop(name)

    def _set_parent(self, parent: Type["_Node"]):
        """
        set the parent and update path
        """

        self._parent = parent
        self._path = f"{self._parent._get_path()}.{self._name}"

    def _get_child(self, name: str):
        """
        return a child object
        """

        return self._children[name]

    def _has_child(self, name: str):
        return name in self._children

    def _get_children(self):
        """
        return a tuple of children
        """

        return tuple(self._children.values())

    def _get_child_from_path(self, path: str):
        """
        get a child from a string path
        """
        nodes = path.split('.')
        _next = self
        for node in nodes:
            _next = _next._get_child(node)

        return _next

    def _get_parent(self):
        return self._parent

    def _get_path(self):

        if self._parent is not None:
            return self._path

        else:
            return self._name

    @property
    def path(self):
        return self._get_path()

    def _update_path(self):

        # recursively update the path

        for name, child in self._children.items():

            child._path = f"{child._parent._get_path()}.{child._name}"
            child._update_path()

    def _change_name(self, name: str):

        self._name = name

        if self._parent is not None:

            self._parent._children.pop(self._name)
            self._set_parent(self._parent)

        # update all the children
        self._update_path()

    @property
    def is_root(self):
        """
        is this the root of the tree
        """
        return self._parent is None

    @property
    def is_leaf(self):
        """
        is this a a leaf of the tree
        """
        if len(self._children) == 0:
            return True

        else:

            return False

    @property
    def name(self):
        return self._name

    def __dir__(self):
        return super().__dir__() + list(self._children.keys())

    def __getattr__(self, name):
        if name in self._children:
            return self._children[name]
        else:
            return super().__getattr__(name)

    def __setattr__(self, name, value):

        if "_children" in self.__dict__:
            if name in self._children:
                raise AttributeError()
            else:
                return super().__setattr__(name, value)
        else:

            return super().__setattr__(name, value)

    def plot_tree(self):

        print(self._plot_tree_console())

    def _plot_tree_console(self):

        if self.is_leaf:
            return self.name

        child_strs = [child._plot_tree_console() for child in list(self._children.values())]
        child_widths = [block_width(s) for s in child_strs]

        # How wide is this block?
        display_width = max(len(self.name),
                            sum(child_widths) + len(child_widths) - 1)

        # Determines midpoints of child blocks
        child_midpoints = []
        child_end = 0
        for width in child_widths:
            child_midpoints.append(child_end + (width // 2))
            child_end += width + 1

        # Builds up the brace, using the child midpoints
        brace_builder = []
        for i in range(display_width):
            if i < child_midpoints[0] or i > child_midpoints[-1]:
                brace_builder.append(' ')
            elif i in child_midpoints:
                brace_builder.append('+')
            else:
                brace_builder.append('-')
        brace = ''.join(brace_builder)

        name_str = '{:^{}}'.format(self.name, display_width)
        below = stack_str_blocks(child_strs)

        return name_str + '\n' + brace + '\n' + below


class Node(_Node):

    # This apparently dumb constructor is needed otherwise pickle will fail

    def __init__(self, name):

        assert is_valid_variable_name(name), "Illegal characters in name %s. You can only use letters and numbers, " \
                                             "and _" % name

        assert name != "name", "You cannot call a node 'name', it is reserved."

        _Node.__init__(self, name)

        #########################################################################

    # This is used by dir() and by the autocompletion in Ipython

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


def block_width(block):
    try:
        return block.index('\n')
    except ValueError:
        return len(block)


def stack_str_blocks(blocks):
    """Takes a list of multiline strings, and stacks them horizontally.

    For example, given 'aaa\naaa' and 'bbbb\nbbbb', it returns
    'aaa bbbb\naaa bbbb'.  As in:

    'aaa  +  'bbbb   =  'aaa bbbb
     aaa'     bbbb'      aaa bbbb'

    Each block must be rectangular (all lines are the same length), but blocks
    can be different sizes.
    """
    builder = []
    block_lens = [block_width(bl) for bl in blocks]
    split_blocks = [bl.split('\n') for bl in blocks]

    for line_list in itertools.zip_longest(*split_blocks, fillvalue=None):
        for i, line in enumerate(line_list):
            if line is None:
                builder.append(' ' * block_lens[i])
            else:
                builder.append(line)
            if i != len(line_list) - 1:
                builder.append(' ')  # Padding
        builder.append('\n')

    return ''.join(builder[:-1])
