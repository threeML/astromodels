from builtins import zip

__author__ = "giacomov"

import collections
import os
import warnings
from dataclasses import dataclass
from typing import Any, Dict, Iterable, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
import scipy.integrate

from astromodels.core.memoization import use_astromodels_memoization
from astromodels.core.my_yaml import my_yaml
from astromodels.core.parameter import IndependentVariable, Parameter
from astromodels.core.property import FunctionProperty
from astromodels.core.tree import DuplicatedNode, Node
from astromodels.functions.function import Function, get_function
from astromodels.sources import (ExtendedSource, ParticleSource, PointSource,
                                 Source, SourceType)
from astromodels.utils.disk_usage import disk_usage
from astromodels.utils.logging import setup_logger
from astromodels.utils.long_path_formatter import long_path_formatter

log = setup_logger(__name__)


class ModelFileExists(IOError):
    pass


class InvalidInput(ValueError):
    pass


class CannotWriteModel(IOError):
    def __init__(self, directory, message):
        # Add a report on disk usage to the message

        free_space = disk_usage(directory).free

        message += "\nFree space on the file system hosting %s was %.2f Mbytes" % (
            directory,
            free_space / 1024.0 / 1024.0,
        )

        super(CannotWriteModel, self).__init__(message)


class ModelInternalError(ValueError):

    pass

@dataclass(frozen=True)
class _LinkedFunctionContainer:
    function_name: str
    linked_path: str
    internal_name: str

    def output(self, rich=True):

        this_dict = collections.OrderedDict()

        tmp = collections.OrderedDict()
        tmp["linked function"] = self.linked_path
        tmp["internal name"] = self.internal_name
        this_dict[self.function_name] = tmp

        output = pd.DataFrame.from_dict(this_dict)

        if rich:

            return output._repr_html_()

        else:

            return output.__repr__()



class Model(Node):
    def __init__(self, *sources):

        # Setup the node, using the special name '__root__' to indicate that this is the root of the tree

        super(Model, self).__init__("__root__")

        # Dictionary to keep point sources

        self._point_sources: Dict[str, PointSource] = collections.OrderedDict()

        # Dictionary to keep extended sources

        self._extended_sources: Dict[
            str, ExtendedSource] = collections.OrderedDict()

        # Dictionary to keep particle sources

        self._particle_sources: Dict[
            str, ParticleSource] = collections.OrderedDict()

        # Loop over the provided sources and process them

        for source in sources:

            self._add_source(source)

        # Now make the list of all the existing parameters

        self._update_parameters()

        # This controls the verbosity of the display
        self._complete_display = False

        # This will keep track of independent variables (if any)
        self._independent_variables = {}

    def _add_source(
            self, source: Union[PointSource, ExtendedSource,
                                ParticleSource]) -> None:
        """
        Remember to call _update_parameters after this!
        :param source:
        :return:
        """

        try:

            self._add_child(source)

        except AttributeError:

            if isinstance(source, Source):

                log.error(
                    "More than one source with the name '%s'. You cannot use the same name for multiple "
                    "sources" % source.name)

                raise DuplicatedNode()

            else:  # pragma: no cover

                raise

        # Now see if this is a point or extended source, and add them to the
        # appropriate dictionary

        if source.source_type == SourceType.POINT_SOURCE:

            self._point_sources[source.name] = source

        elif source.source_type == SourceType.EXTENDED_SOURCE:

            self._extended_sources[source.name] = source

        elif source.source_type == SourceType.PARTICLE_SOURCE:

            self._particle_sources[source.name] = source

        else:  # pragma: no cover

            raise InvalidInput(
                "Input sources must be either a point source or an extended source"
            )

    def _remove_source(self, source_name: str):
        """
        Remember to call _update_parameters after this
        :param source_name:
        :return:
        """

        if not source_name in self.sources:

            log.error(f"Source {source_name} is not part of the current model")

            raise AssertionError()

        source = self.sources.pop(source_name)

        if source.source_type == SourceType.POINT_SOURCE:

            self._point_sources.pop(source.name)

        elif source.source_type == SourceType.EXTENDED_SOURCE:

            self._extended_sources.pop(source.name)

        elif source.source_type == SourceType.PARTICLE_SOURCE:

            self._particle_sources.pop(source.name)

        self._remove_child(source_name)

    def _find_parameters(self, node) -> Dict[str, Parameter]:
        
        return self._recursively_gather_node_type(self, Parameter)

    def _find_properties(self, node) -> Dict[str, FunctionProperty]:

        return self._recursively_gather_node_type(self, FunctionProperty)

    
    def _update_parameters(self) -> None:

        self._parameters: Dict[str, Parameter] = self._find_parameters(self)
        self._properties: Dict[str, FunctionProperty] = self._find_properties(self)
        
    @property
    def parameters(self) -> Dict[str, Parameter]:
        """
        Return a dictionary with all parameters

        :return: dictionary of parameters
        """
        self._update_parameters()

        return self._parameters

    @property
    def free_parameters(self) -> Dict[str, Parameter]:
        """
        Get a dictionary with all the free parameters in this model

        :return: dictionary of free parameters
        """

        # Refresh the list

        self._update_parameters()

        # Filter selecting only free parameters

        free_parameters_dictionary: Dict[
            str, Parameter] = collections.OrderedDict()

        for parameter_name, parameter in list(self._parameters.items()):

            if parameter.free:

                free_parameters_dictionary[parameter_name] = parameter

        return free_parameters_dictionary

    @property
    def linked_parameters(self) -> Dict[str, Parameter]:
        """
        Get a dictionary with all parameters in this model in a linked status. A parameter is in a linked status
        if it is linked to another parameter (i.e. it is forced to have the same value of the other parameter), or
        if it is linked with another parameter or an independent variable through a law.

        :return: dictionary of linked parameters
        """

        # Refresh the list

        self._update_parameters()

        # Filter selecting only free parameters

        linked_parameter_dictionary: Dict[
            str, Parameter] = collections.OrderedDict()

        for parameter_name, parameter in list(self._parameters.items()):

             if parameter.has_auxiliary_variable:

                linked_parameter_dictionary[parameter_name] = parameter

        return linked_parameter_dictionary

    @property
    def properties(self) -> Dict[str, Parameter]:
        """
        Return a dictionary with all parameters

        :return: dictionary of parameters
        """
        self._update_parameters()

        return self._properties

    
    @property
    def linked_functions(self) -> List[_LinkedFunctionContainer]:
        """
        return a list of containers for the linked functions
        
        """

        linked_functions = []

        for function in self._get_all_functions():

            
            
            if "composite" in function._children:

                data = function.to_dict()['composite']
                
                # this is a composite

                if "external_functions" in data:

                    for k, v in data["external_functions"].items():

                        if v:

                            for name, path in v.items():

                                tmp = _LinkedFunctionContainer(function_name=function.path,
                                                               linked_path=path,
                                                               internal_name=name)

                                linked_functions.append(tmp)
            else:

                data = function.to_dict()[function.shape.name]
                
                if "external_functions" in data:

                    for name, path in data["external_functions"].items():
                        
                        tmp = _LinkedFunctionContainer(function_name=function.path,
                                                       linked_path=path,
                                                       internal_name=name)

                
                        linked_functions.append(tmp)


        return linked_functions

    def _get_all_functions(self) -> List[Function]:
        all_functions = []

        for source_name, source in self.sources.items():

            # look at the spectrum node

            all_functions.extend(source.spectrum._get_children())

        return all_functions

    
    def set_free_parameters(self, values: Iterable[float]) -> None:
        """
        Set the free parameters in the model to the provided values.

        NOTE: of course, order matters

        :param values: a list of new values
        :return: None
        """

        if not len(values) == len(self.free_parameters):

            log.error(
                f"tried to pass {len(values)} parameters but need {len(self.free_parameters)}"
            )

            raise AssertionError()

        for parameter, this_value in zip(list(self.free_parameters.values()),
                                         values):

            parameter.value = this_value

    def __getitem__(self, path: str):
        """
        Get a parameter from a path like "source_1.component.powerlaw.logK". This might be useful in certain
        context, although in an interactive analysis there is no reason to use this.

        :param path: the address of the parameter
        :return: the parameter
        """

        return self._get_child_from_path(path)

    def __contains__(self, path: str) -> bool:
        """
        This allows the model to be used with the "in" operator, like;

        > if 'myparameter' in model:
        >    print("Myparameter is contained in the model")

        :param path: the parameter to look for
        :return:
        """

        try:

            _ = self._get_child_from_path(path)

        except (AttributeError, KeyError, TypeError):

            return False

        else:

            return True

    def __iter__(self):
        """
        This allows the model to be iterated on, like in:

        for parameter in model:
            ...

        NOTE: this will iterate over *all* parameters in the model, also those that are not free (and thus are not
        normally displayed). If you need to operate only on free parameters, just check if they are free within
        the loop or use the .free_parameters dictionary directly

        :return: iterator
        """

        for parameter in self.parameters:

            yield self.parameters[parameter]

    @property
    def point_sources(self) -> Dict[str, PointSource]:
        """
        Returns the dictionary of all defined point sources

        :return: collections.OrderedDict()
        """
        return self._point_sources

    @property
    def extended_sources(self) -> Dict[str, ExtendedSource]:
        """
        Returns the dictionary of all defined extended sources

        :return: collections.OrderedDict()

        """
        return self._extended_sources

    @property
    def particle_sources(self) -> Dict[str, ParticleSource]:
        """
        Returns the dictionary of all defined particle sources

        :return: collections.OrderedDict()

        """
        return self._particle_sources

    @property
    def sources(
            self
    ) -> Dict[str, Union[PointSource, ExtendedSource, ParticleSource]]:
        """
        Returns a dictionary containing all defined sources (of any kind)

        :return: collections.OrderedDict()

        """

        sources: Dict[str, Union[PointSource, ExtendedSource,
                                 ParticleSource]] = collections.OrderedDict()

        for d in (self.point_sources, self.extended_sources,
                  self.particle_sources):

            sources.update(d)

        return sources

    def add_source(
        self, new_source: Union[PointSource, ExtendedSource,
                                 ParticleSource]) -> None:
        """
        Add the provided source to the model

        :param new_source: the new source to be added (an instance of PointSource, ExtendedSource or ParticleSource)
        :return: (none)
        """

        self._add_source(new_source)

        self._update_parameters()

    def remove_source(self, source_name: str) -> None:
        """
        Returns a new model with the provided source removed from the current model. Any parameters linked to the source to be removed are automatically unlinked.

        :param source_name: the name of the source to be removed
        :return: a new Model instance without the source
        """

        self.unlink_all_from_source(source_name, warn=True)

        self._remove_source(source_name)

        self._update_parameters()

    def unlink_all_from_source(self,
                               source_name: str,
                               warn: bool = False) -> None:
        """
        Unlink all parameters of the current model that are linked to a parameter of a given source.
        To be called before removing a source from the model.

        :param source_name: the name of the source to which to remove all links
        :param warn: If True, prints a warning if any parameters were unlinked.
        """

        tempmodel = Model(self[source_name])
        unlinked_parameters = collections.OrderedDict()

        for par in self.linked_parameters.values():

            target = par._aux_variable['variable']

            if target.path in tempmodel:

                unlinked_parameters[par.name] = par
                self.unlink(par)

        if warn and unlinked_parameters:

            log.warning(
                "The following %d parameters that were linked to source %s have been automatically un-linked: %s"
                % (len(unlinked_parameters), source_name,
                   [p.path for p in unlinked_parameters.values()]),
                RuntimeWarning)

    def add_independent_variable(self, variable: IndependentVariable) -> None:
        """
        Add a global independent variable to this model, such as time.

        :param variable: an IndependentVariable instance
        :return: none
        """

        if not isinstance(variable, IndependentVariable):
            log.error("Variable must be an instance of IndependentVariable")

            raise AssertionError()

        if self._has_child(variable.name):

            self._remove_child(variable.name)

        self._add_child(variable)

        # Add also to the list of independent variables
        self._independent_variables[variable.name] = variable

    def remove_independent_variable(self, variable_name: str) -> None:
        """
        Remove an independent variable which was added with add_independent_variable

        :param variable_name: name of variable to remove
        :return:
        """

        self._remove_child(variable_name)

        # Remove also from the list of independent variables
        self._independent_variables.pop(variable_name)

    def add_external_parameter(self, parameter: Parameter) -> None:
        """
        Add a parameter that comes from something other than a function, to the model.

        :param parameter: a Parameter instance
        :return: none
        """

        if not isinstance(parameter, Parameter):

            log.error("Variable must be an instance of IndependentVariable")

            raise AssertionError()

        if self._has_child(parameter.name):

            # Remove it from the children only if it is a Parameter instance, otherwise don't, which will
            # make the _add_child call fail (which is the expected behaviour! You shouldn't call two children
            # with the same name)

            if isinstance(self._get_child(parameter.name), Parameter):

                log.warning(
                    "External parameter %s already exist in the model. Overwriting it..."
                    % parameter.name)

                self._remove_child(parameter.name)

        # This will fail if another node with the same name is already in the model

        self._add_child(parameter)

    def remove_external_parameter(self, parameter_name: str) -> None:
        """
        Remove an external parameter which was added with add_external_parameter

        :param variable_name: name of parameter to remove
        :return:
        """

        self._remove_child(parameter_name)

    def link(self, parameter_1, parameter_2, link_function=None) -> None:
        """
        Link the value of the provided parameters through the provided function (identity is the default, i.e.,
        parameter_1 = parameter_2).

        :param parameter_1: the first parameter;can be either a single parameter or a list of prarameters
        :param parameter_2: the second parameter
        :param link_function: a function instance. If not provided, the identity function will be used by default.
        Otherwise, this link will be set: parameter_1 = link_function(parameter_2)
        :return: (none)
        """
        if not isinstance(parameter_1, list):
            # Make a list of one element
            parameter_1_list = [parameter_1]
        else:
            # Make a copy to avoid tampering with the input
            parameter_1_list = list(parameter_1)

        for param_1 in parameter_1_list:
            if not param_1.path in self:

                log.error( "Parameter %s is not contained in this model" % param_1.path )

                raise AssertionError()
                
        if not parameter_2.path in self:

            log.error( "Parameter %s is not contained in this model" % parameter_2.path )

            raise AssertionError()
            
        if link_function is None:
            # Use the Line function by default, with both parameters fixed so that the two
            # parameters to be linked will vary together
            link_function = get_function("Line")

            link_function.a.value = 0
            link_function.a.fix = True

            link_function.b.value = 1
            link_function.b.fix = True

        for param_1 in parameter_1_list:
            param_1.add_auxiliary_variable(parameter_2, link_function)
            # Now set the units of the link function
            link_function.set_units(parameter_2.unit, param_1.unit)

    def unlink(self, parameter: Parameter) -> None:
        """
        Sets free one or more parameters which have been linked previously

        :param parameter: the parameter to be set free, can also be a list of parameters
        :return: (none)
        """

        if not isinstance(parameter, list):
            # Make a list of one element
            parameter_list = [parameter]
        else:
            # Make a copy to avoid tampering with the input
            parameter_list = list(parameter)

        for param in parameter_list:
            if param.has_auxiliary_variable:
                param.remove_auxiliary_variable()

            else:

                with warnings.catch_warnings():

                    warnings.simplefilter("always", RuntimeWarning)

                    log.warning(
                        "Parameter %s has no link to be removed." % param.path
                       
                    )

    def display(self, complete: bool=False) -> None:
        """
        Display information about the point source.

        :param complete : if True, displays also information on fixed parameters
        :return: (none)
        """

        # Switch on the complete display flag
        self._complete_display = bool(complete)

        # This will automatically choose the best representation among repr and repr_html

        super(Model, self).display()

        # Go back to default

        self._complete_display = False

    def _repr__base(self, rich_output=False):

        if rich_output:

            new_line = "<br>"

        else:

            new_line = "\n"

        # Table with the summary of the various kind of sources
        sources_summary = pd.DataFrame.from_dict(
            collections.OrderedDict(
                [
                    ("Point sources", [self.get_number_of_point_sources()]),
                    ("Extended sources", [self.get_number_of_extended_sources()]),
                    ("Particle sources", [self.get_number_of_particle_sources()]),
                ]
            ),
            columns=["N"],
            orient="index",
        )

        # These properties traverse the whole tree everytime, so let's cache their results here
        parameters = self.parameters
        free_parameters = self.free_parameters
        linked_parameters = self.linked_parameters
        properties = self.properties
        linked_functions = self.linked_functions
        
        # Summary of free parameters
        if len(free_parameters) > 0:

            parameter_dict = collections.OrderedDict()

            for parameter_name, parameter in list(free_parameters.items()):
                # Generate table with only a minimal set of info

                # Generate table with only a minimal set of info
                if rich_output:

                    this_name = long_path_formatter(parameter_name, 70)

                else:

                    # In a terminal we need to use less characters

                    this_name = long_path_formatter(parameter_name, 40)

                d = parameter.to_dict()
                parameter_dict[this_name] = collections.OrderedDict()

                for key in ["value", "unit", "min_value", "max_value"]:

                    parameter_dict[this_name][key] = d[key]

            free_parameters_summary = pd.DataFrame.from_dict(parameter_dict).T

            # Re-order it
            free_parameters_summary = free_parameters_summary[
                ["value", "min_value", "max_value", "unit"]
            ]

        else:

            free_parameters_summary = pd.DataFrame()

        if len(parameters) - len(free_parameters) - len(linked_parameters) > 0:

            fixed_parameter_dict = collections.OrderedDict()

            for parameter_name, parameter in list(parameters.items()):

                if parameter.free or parameter_name in linked_parameters:

                    continue

                # Generate table with only a minimal set of info
                if rich_output:

                    this_name = long_path_formatter(parameter_name, 70)

                else:

                    # In a terminal we need to use less characters

                    this_name = long_path_formatter(parameter_name, 40)

                d = parameter.to_dict()
                fixed_parameter_dict[this_name] = collections.OrderedDict()

                for key in ["value", "unit", "min_value", "max_value"]:

                    fixed_parameter_dict[this_name][key] = d[key]

            fixed_parameters_summary = pd.DataFrame.from_dict(fixed_parameter_dict).T

            # Re-order it
            fixed_parameters_summary = fixed_parameters_summary[
                ["value", "min_value", "max_value", "unit"]
            ]

        else:

            fixed_parameters_summary = pd.DataFrame()

        # Summary of linked parameters

        linked_frames = []

        if linked_parameters:

            for parameter_name, parameter in list(linked_parameters.items()):

                parameter_dict = collections.OrderedDict()

                # Generate table with only a minimal set of info

                variable, law = parameter.auxiliary_variable

                this_dict = collections.OrderedDict()

                this_dict["linked to"] = variable.path
                this_dict["function"] = law.name
                this_dict["current value"] = parameter.value
                this_dict["unit"] = parameter.unit

                parameter_dict[parameter_name] = this_dict

                this_parameter_frame = pd.DataFrame.from_dict(parameter_dict)

                linked_frames.append(this_parameter_frame)

        else:

            # No linked parameters

            pass

        
        
        empty_frame = "(none)%s" % new_line

        # Summary of free parameters
        if len(properties) > 0:

            property_dict = collections.OrderedDict()

            for property_name, prop in properties.items():

                # Generate table with only a minimal set of info
                if rich_output:

                    this_name = long_path_formatter(property_name, 70)

                else:

                    # In a terminal we need to use less characters

                    this_name = long_path_formatter(property_name, 40)

                d = prop.to_dict()
                property_dict[this_name] = collections.OrderedDict()
                
                for key in ["value", "allowed values"]:

                    property_dict[this_name][key] = d[key]

            properties_summary = pd.DataFrame.from_dict(property_dict).T

            # Re-order it
            properties_summary = properties_summary[
                ["value", "allowed values"]
            ]

        else:

            properties_summary = pd.DataFrame()

        empty_frame = "(none)%s" % new_line
        
        # Independent variables

        independent_v_frames = []

        if self._independent_variables:

            for variable_name, variable_instance in list(
                self._independent_variables.items()
            ):

                v_dict = collections.OrderedDict()

                # Generate table with only a minimal set of info

                this_dict = collections.OrderedDict()

                this_dict["current value"] = variable_instance.value
                this_dict["unit"] = variable_instance.unit

                v_dict[variable_name] = this_dict

                this_v_frame = pd.DataFrame.from_dict(v_dict)

                independent_v_frames.append(this_v_frame)

        else:

            # No independent variables

            pass


        
        if rich_output:

            source_summary_representation = sources_summary._repr_html_()

            if free_parameters_summary.empty:

                free_parameters_representation = empty_frame

            else:

                free_parameters_representation = free_parameters_summary._repr_html_()

            if len(linked_frames) == 0:

                linked_summary_representation = empty_frame

            else:

                linked_summary_representation = ""

                for linked_frame in linked_frames:

                    linked_summary_representation += linked_frame._repr_html_()
                    linked_summary_representation += new_line

            if properties_summary.empty:

                properties_representation = empty_frame

            else:

                properties_representation = properties_summary._repr_html_()

                    
            if len(linked_functions) == 0:

                linked_function_summary_representation = empty_frame

            else:

                linked_function_summary_representation = ""

                for linked_function in linked_functions:

                    linked_function_summary_representation += linked_function.output(rich=True)
                    linked_function_summary_representation += new_line


                    
            if len(independent_v_frames) == 0:

                independent_v_representation = empty_frame

            else:

                independent_v_representation = ""

                for v_frame in independent_v_frames:

                    independent_v_representation += v_frame._repr_html_()
                    independent_v_representation += new_line

            if fixed_parameters_summary.empty:

                fixed_parameters_representation = empty_frame

            else:

                fixed_parameters_representation = fixed_parameters_summary._repr_html_()

        else:

            source_summary_representation = sources_summary.__repr__()

            if free_parameters_summary.empty:

                free_parameters_representation = empty_frame

            else:

                free_parameters_representation = free_parameters_summary.__repr__()

            if properties_summary.empty:

                properties_representation = empty_frame

            else:

                properties_representation = properties_summary.__repr__()


                
            if len(linked_frames) == 0:

                linked_summary_representation = empty_frame

            else:

                linked_summary_representation = ""

                for linked_frame in linked_frames:

                    linked_summary_representation += linked_frame.__repr__()
                    linked_summary_representation += "%s%s" % (new_line, new_line)

            if len(linked_functions) == 0:

                linked_function_summary_representation = empty_frame

            else:

                linked_function_summary_representation = ""

                for linked_function in linked_functions:

                    linked_function_summary_representation += linked_function.output(rich=False)
                    linked_function_summary_representation += "%s%s" % (new_line, new_line)

                    
            if len(independent_v_frames) == 0:

                independent_v_representation = empty_frame

            else:

                independent_v_representation = ""

                for v_frame in independent_v_frames:

                    independent_v_representation += v_frame.__repr__()
                    independent_v_representation += "%s%s" % (new_line, new_line)

            if fixed_parameters_summary.empty:

                fixed_parameters_representation = empty_frame

            else:

                fixed_parameters_representation = fixed_parameters_summary.__repr__()

        # Build the representation

        representation = "Model summary:%s" % (new_line)

        if not rich_output:

            representation += "==============%s%s" % (new_line, new_line)

        else:

            representation += new_line

        # Summary on sources

        representation += source_summary_representation

        representation += new_line

        # Free parameters

        representation += "%sFree parameters (%i):%s" % (
            new_line,
            len(free_parameters),
            new_line,
        )

        if not rich_output:

            representation += "--------------------%s%s" % (new_line, new_line)

        else:

            representation += new_line

        representation += free_parameters_representation

        representation += new_line

        # Fixed parameters

        n_fix = len(parameters) - len(free_parameters) - len(linked_parameters)

        representation += "%sFixed parameters (%i):%s" % (new_line, n_fix, new_line)

        if self._complete_display:

            if not rich_output:

                representation += "---------------------%s%s" % (new_line, new_line)

            else:

                representation += new_line

            representation += fixed_parameters_representation

        else:

            representation += (
                "(abridged. Use complete=True to see all fixed parameters)%s" % new_line
            )

        representation += new_line

       # Properties

        representation += "%sProperties (%i):%s" % (
            new_line,
            len(properties),
            new_line,
        )

        if not rich_output:

            representation += "--------------------%s%s" % (new_line, new_line)

        else:

            representation += new_line

        representation += properties_representation

        representation += new_line

        
        # Linked parameters

        representation += "%sLinked parameters (%i):%s" % (
            new_line,
            len(self.linked_parameters),
            new_line,
        )

        if not rich_output:

            representation += "----------------------%s%s" % (new_line, new_line)

        else:

            representation += new_line

        representation += linked_summary_representation

        # Independent variables

        representation += "%sIndependent variables:%s" % (new_line, new_line)

        if not rich_output:

            representation += "----------------------%s%s" % (new_line, new_line)

        else:

            representation += new_line

        representation += independent_v_representation

        # Linked functions

        representation += "%sLinked functions (%i):%s" % (
            new_line,
            len(self.linked_functions),
            new_line,
        )

        if not rich_output:

            representation += "----------------------%s%s" % (new_line, new_line)

        else:

            representation += new_line

        representation += linked_function_summary_representation

        
        return representation

    def to_dict_with_types(self):

        # Get the serialization dictionary

        data = self.to_dict()

        # Add the types to the sources

        for key in list(data.keys()):

            try:

                element = self._get_child(key)

            except KeyError:  # pragma: no cover

                raise RuntimeError("Source %s is unknown" % key)

            else:

                # There are three possible cases. Either the element is a source, or it is an independent
                # variable, or a parameter

                if hasattr(element, "source_type"):

                    # Change the name of the key adding the source type

                    data["%s (%s)" % (key, element.source_type)] = data.pop(key)

                elif isinstance(element, IndependentVariable):

                    data["%s (%s)" % (key, "IndependentVariable")] = data.pop(key)

                elif isinstance(element, Parameter):

                    data["%s (%s)" % (key, "Parameter")] = data.pop(key)

                # elif isinstance(element, FunctionProperty):

                #     data["%s (%s)" % (key, "FunctionProperty")] = data.pop(key)

                    
                else:  # pragma: no cover

                    raise ModelInternalError("Found an unknown class at the top level")

        return data

    def save(self, output_file, overwrite=False):

        """Save the model to disk"""

        if os.path.exists(output_file) and overwrite is False:

            raise ModelFileExists(
                "The file %s exists already. If you want to overwrite it, use the 'overwrite=True' "
                "options as 'model.save(\"%s\", overwrite=True)'. "
                % (output_file, output_file)
            )

        else:

            data = self.to_dict_with_types()

            # Write it to disk

            try:

                # Get the YAML representation of the data

                representation = my_yaml.dump(data, default_flow_style=False)

                with open(output_file, "w+") as f:

                    # Add a new line at the end of each voice (just for clarity)

                    f.write(representation.replace("\n", "\n\n"))

            except IOError:

                raise CannotWriteModel(
                    os.path.dirname(os.path.abspath(output_file)),
                    "Could not write model file %s. Check your permissions to write or the "
                    "report on the free space which follows: " % output_file,
                )

    def get_number_of_point_sources(self) ->int:
        """
        Return the number of point sources

        :return: number of point sources
        """
        return len(self._point_sources)

    def get_point_source_position(self, id) -> Tuple[float]:
        """
        Get the point source position (R.A., Dec)

        :param id: id of the source
        :return: a tuple with R.A. and Dec.
        """

        pts = list(self._point_sources.values())[id]

        return pts.position.get_ra(), pts.position.get_dec()

    def get_point_source_fluxes(self, id: int, energies: np.ndarray, tag=None) -> np.ndarray:
        """
        Get the fluxes from the id-th point source

        :param id: id of the source
        :param energies: energies at which you need the flux
        :param tag: a tuple (integration variable, a, b) specifying the integration to perform. If this
        parameter is specified then the returned value will be the average flux for the source computed as the integral
        between a and b over the integration variable divided by (b-a). The integration variable must be an independent
        variable contained in the model. If b is None, then instead of integrating the integration variable will be
        set to a and the model evaluated in a.
        :return: fluxes
        """
        
        return list(self._point_sources.values())[id](energies, tag=tag)

    def get_point_source_name(self, id: int) -> str:

        return list(self._point_sources.values())[id].name

    def get_number_of_extended_sources(self) -> int:
        """
        Return the number of extended sources

        :return: number of extended sources
        """
        return len(self._extended_sources)

    def get_extended_source_fluxes(self, id: int, j2000_ra: float, j2000_dec: float, energies: np.ndarray) -> np.ndarray:
        """
        Get the flux of the id-th extended sources at the given position at the given energies

        :param id: id of the source
        :param j2000_ra: R.A. where the flux is desired
        :param j2000_dec: Dec. where the flux is desired
        :param energies: energies at which the flux is desired
        :return: flux array
        """

        return list(self._extended_sources.values())[id](j2000_ra, j2000_dec, energies)

    def get_extended_source_name(self, id: int) -> str:
        """
        Return the name of the n-th extended source

        :param id: id of the source (integer)
        :return: the name of the id-th source
        """

        return list(self._extended_sources.values())[id].name

    def get_extended_source_boundaries(self, id: int):

        (ra_min, ra_max), (dec_min, dec_max) = list(self._extended_sources.values())[
            id
        ].get_boundaries()

        return ra_min, ra_max, dec_min, dec_max

    def is_inside_any_extended_source(self, j2000_ra: float, j2000_dec: float) -> bool:

        for ext_source in list(self.extended_sources.values()):

            (ra_min, ra_max), (dec_min, dec_max) = ext_source.get_boundaries()

            # Transform from the 0...360 convention to the -180..180 convention, so that
            # the comparison is easier
            if ra_min > 180:

                ra_min = -(360 - ra_min)

            if ra_min <= j2000_ra <= ra_max and dec_min <= j2000_dec <= dec_max:

                return True

        # If we are here, it means that no extended source contains the provided coordinates

        return False

    def get_number_of_particle_sources(self) -> int:
        """
        Return the number of particle sources

        :return: number of particle sources
        """

        return len(self._particle_sources)

    def get_particle_source_fluxes(self, id: int, energies: np.ndarray) -> np.ndarray:
        """
        Get the fluxes from the id-th point source

        :param id: id of the source
        :param energies: energies at which you need the flux
        :return: fluxes
        """

        return list(self._particle_sources.values())[id](energies)

    def get_particle_source_name(self, id: int) -> str:

        return list(self._particle_sources.values())[id].name

    def get_total_flux(self, energies: np.ndarray) -> float:
        """
        Returns the total differential flux at the provided energies from all *point* sources

        :return:
        """

        fluxes = []

        for src in self._point_sources:

            fluxes.append(self._point_sources[src](energies))

        return np.sum(fluxes, axis=0)

