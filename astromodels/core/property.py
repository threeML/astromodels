import collections
import copy
from typing import Any, Dict, List, Optional, Tuple, Union

from astromodels.utils.logging import setup_logger

from .tree import Node

log = setup_logger(__name__)

# Exception for when a parameter is out of its bounds
class SettingUnknownValue(RuntimeError):
    pass


class PropertyBase(Node):

    def __init__(self,
                 name: str,
                 desc: str,
                 value: Optional[str] = None,
                 allowed_values: Optional[List[str]] = None,
                 defer: bool = False,
                 eval_func: Optional[str] = None
                 ):

        # Make this a node

        Node.__init__(self, name)

        self._allowed_values: Optional[List[str]] = allowed_values
        self._defer: bool = defer
        self._eval_func: Optional[str] = eval_func
        
        if (value is None) and (not self._defer):

            log.error(f"property {name} was given no initial value but is NOT deferred")

        # now we set the value

        self.value = value

        self.__doc__ = desc

        self._desc = desc

    def _get_value(self) -> Any:
        """
        Return current parameter value
        """

        log.debug(f"accessing the property {self.name} with value {self._internal_value}")
        
        return self._internal_value
    

    def _set_value(self, new_value) -> None:
        """
        Sets the current value of the parameter, ensuring that it is within the allowed range.
        """

        if (self._defer) and (new_value is None):

                # this is ok
                pass

        elif self._allowed_values is not None:

            if not new_value in self._allowed_values:

                log.error(f"{self.name} can only take the values {','.join(self._allowed_values)} not {new_value}")

                raise SettingUnknownValue()
                

        self._internal_value = new_value
                                       
        # if there is an eval func value
        # then we need to execute the function
        # on the parent

        if (self._internal_value == "_tmp") and self._defer:

            # do not execute in this mode
            
            return 

        if self._eval_func is not None:

            # if there is a parent
            if self._parent is not None:

                if self._parent.name == "composite":
                    # ok, we have a composite function

                    func_idx = int(self._name.split("_")[-1]) - 1

                    getattr(self._parent._functions[func_idx], str(self._eval_func))()

                else:
                
                    getattr(self._parent, str(self._eval_func))()

            # other wise this will run when the parent is set
                
        
        

    value = property(
        _get_value,
        _set_value,
        doc=
        "Get and sets the current value for the propert",)

    def _set_parent(self, parent):

        # we intecept here becuase we want
        # to make sure the eval works
        
        super(PropertyBase, self)._set_parent(parent)

        # now we want to update because we have a parent
        self.value = self._internal_value
        

    
    @property
    def is_deferred(self) -> bool:
        return self._defer

    @property
    def description(self) -> Optional[str]:
        """
        Return a description of this parameter

        :return: a string cointaining a description of the meaning of this parameter
        """
        return self._desc
    
    def duplicate(self) -> "FunctionProperty":
        """
        Returns an exact copy of the current property
        """

        # Deep copy everything to make sure that there are no ties between the new instance and the old one

        new_property = copy.deepcopy(self)

        return new_property

    def _repr__base(self, rich_output):  # pragma: no cover

        raise NotImplementedError(
            "You need to implement this for the actual Property class")


    @staticmethod
    def _to_python_type(variable):
        """
        Returns the value in the variable handling also np.array of one element

        :param variable: input variable
        :return: the value of the variable having a python type (int, float, ...)
        """

        # Assume variable is a np.array, fall back to the case where variable is already a primitive type

        try:

            return variable.item()

        except AttributeError:

            return variable

    
    def to_dict(self, minimal=False) -> Dict[str, Any]:
        """Returns the representation for serialization"""

        data = collections.OrderedDict()

        if minimal:

            # In the minimal representation we just output the value

            data["value"] = self._to_python_type(self.value)

        else:

            # In the complete representation we output everything is needed to re-build the object

            data["value"] = str(self.value)
            data["desc"] = str(self._desc)
            data["allowed values"] = self._to_python_type(self._allowed_values)
            data["defer"] = self._to_python_type(self._defer)
            data["function"] = str(self._eval_func)

        return data

    
        
class FunctionProperty(PropertyBase):

        def __init__(self,
                     name: str,
                     desc: str,
                     value: Optional[str] = None,
                     allowed_values: Optional[List[Any]] = None,
                     defer: bool = False,
                     eval_func: Optional[str] = None
                     
                     ):

            super(FunctionProperty, self).__init__(name=name,desc=desc,
                                                   value=value,
                                                   allowed_values=allowed_values,
                                                   defer=defer,
                                                   eval_func=eval_func)

            
        def _repr__base(self, rich_output=False):

            representation = (
                f"Property {self.name} = {self.value}\n"
                f"(allowed values = {'all' if self._allowed_values is None else ' ,'.join(self._allowed_values)})")
                
            return representation
            
