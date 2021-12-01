import collections
import os
import uuid
from typing import Any, Dict, List, Optional, Union

from future import standard_library

from astromodels.core.node_type import NodeBase
from astromodels.utils.io import display
from astromodels.utils.logging import setup_logger
from astromodels.utils.valid_variable import is_valid_variable_name

standard_library.install_aliases()


class DuplicatedNode(Exception):
    pass


class ProtectedAttribute(RuntimeError):
    pass


class NonExistingAttribute(RuntimeWarning):
    pass


log = setup_logger(__name__)

class Node(NodeBase):

    # This apparently dumb constructor is needed otherwise pickle will fail

    def __init__(self, name: str):

        if not is_valid_variable_name(name):

            log.error("Illegal characters in name %s. You can only use letters and numbers, " \
                                             "and _" % name)

            raise AssertionError()

        if not name != "name":

            log.error("You cannot call a node 'name', it is reserved.")

            raise AssertionError()

        NodeBase.__init__(self, name)
        
        #########################################################################

    def __hash__(self):

        # return the hash of the path
        # anytime we need a hash should
        # be after a model is finalized
        # so it should be static
        
        return hash(self._get_path())
        
        
    def to_dict(self, minimal: bool=False) -> Dict[str, Any]:
        """

        """    
        this_dict: Dict[str, Any] = collections.OrderedDict()

        for child in self._get_children():

            this_dict[child.name] = child.to_dict(minimal)

        return this_dict

    # This is used by dir() and by the autocompletion in Ipython

    def __dir__(self):

        # Get the names of the attributes of the class
        l = list(self.__class__.__dict__.keys())

        # Get all the children
        l.extend([child.name for child in self._get_children()])

        return l

    def _repr__base(self, rich_output):

        return f"Node with name {self.name}"

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
