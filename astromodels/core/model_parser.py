from builtins import object, str

__author__ = "giacomov"

import re
import warnings
from typing import Any, Dict, List, Optional, Union

from astromodels.core import (model, parameter, polarization, sky_direction,
                              spectral_component)
from astromodels.core.my_yaml import my_yaml
from astromodels.functions import function
from astromodels.sources import extended_source, particle_source, point_source
from astromodels.sources.source import SourceType
from astromodels.utils.logging import setup_logger

log = setup_logger(__name__)


class ModelIOError(IOError):
    pass


class ModelYAMLError(my_yaml.YAMLError):
    pass


class ModelSyntaxError(RuntimeError):
    pass


def load_model(filename):
    """
    Load a model from a file.

    :param filename: the name of the file containing the model
    :return: an instance of a Model
    """

    parser = ModelParser(filename)

    return parser.get_model()


def clone_model(model_instance):
    """
    Returns a copy of the given model with all objects cloned. This is equivalent to saving the model to
    a file and reload it, but it doesn't require writing or reading to/from disk. The original model is not touched.

    :param model: model to be cloned
    :return: a cloned copy of the given model
    """

    data = model_instance.to_dict_with_types()

    parser = ModelParser(model_dict=data)

    return parser.get_model()


def model_unpickler(state):

    return ModelParser(model_dict=state).get_model()


class ModelParser(object):
    def __init__(self, model_file=None, model_dict=None):

        if not ((model_file is not None) or (model_dict is not None)):

            log.error("You have to provide either a model file or a "
                      "model dictionary")

            raise AssertionError()

        if model_file is not None:

            # Read model file and deserialize into a dictionary

            try:

                with open(model_file) as f:

                    self._model_dict = my_yaml.load(f,
                                                    Loader=my_yaml.FullLoader)

            except IOError:

                log.error(
                    "File %s cannot be read. Check path and permissions for current user."
                    % model_file)

                raise ModelIOError()

            except my_yaml.YAMLError:

                log.error("Could not parse file %s. Check your syntax." %
                          model_file)

                raise ModelYAMLError()

        else:

            self._model_dict = model_dict

        self._parse()

    def _parse(self):

        # Traverse the dictionary and create all the needed classes

        # The first level is the source level

        self._sources = []
        self._independent_variables = []
        self._external_parameters = []
        self._links = []
        self._external_parameter_links = []
        self._extra_setups = []
        self._external_functions = []
        for source_or_var_name, source_or_var_definition in list(
                self._model_dict.items()):


            # first look for independent variable
            
            if source_or_var_name.find("(IndependentVariable)") > 0:

                var_name = source_or_var_name.split("(")[0].replace(" ", "")

                this_parser = IndependentVariableParser(
                    var_name, source_or_var_definition)

                res = this_parser.get_variable()

                assert isinstance(res, parameter.IndependentVariable)

                self._independent_variables.append(res)

            elif source_or_var_name.find("(Parameter)") > 0:

                var_name = source_or_var_name.split("(")[0].replace(" ", "")

                this_parser = ParameterParser(var_name,
                                              source_or_var_definition)

                res = this_parser.get_variable()

                assert isinstance(res, parameter.Parameter)

                self._external_parameters.append(res)

                self._links.extend(this_parser.links)
            #                self._external_parameter_links.extend(this_parser.links)

            else:

                this_parser = SourceParser(source_or_var_name,
                                           source_or_var_definition)

                res = this_parser.get_source()

                assert (isinstance(res, point_source.PointSource)
                        or isinstance(res, extended_source.ExtendedSource)
                        or isinstance(res, particle_source.ParticleSource))

                self._sources.append(res)

                self._links.extend(this_parser.links)

                self._extra_setups.extend(this_parser.extra_setups)
                
                self._external_functions.extend(this_parser.external_functions)

                
    def get_model(self):

        # Instance the model with all the parsed sources

        new_model = model.Model(*self._sources)

        # Now set up IndependentVariable instances (if any)

        for independent_variable in self._independent_variables:

            new_model.add_independent_variable(independent_variable)

        # Now set up external parameters (if any)
        for parameter in self._external_parameters:

            new_model.add_external_parameter(parameter)

        # Now set up the links

        for link in self._links:

            path = link["parameter_path"]
            variable = link["variable"]
            law = link["law"]

            new_model[path].add_auxiliary_variable(new_model[variable], law)

        # the extra_setups (if any)

        for extra_setup in self._extra_setups:

            path = extra_setup["function_path"]

            for property, value in list(extra_setup["extra_setup"].items()):

                log.debug(f"adding {property} with {value}")

                # First, check to see if the we have a valid path in the new model.
                # If we aren't given a path, interpret it as being given a value.
                if value in new_model:
                    new_model[path].__setattr__(property, new_model[value])
                else:
                    new_model[path].__setattr__(property, value)

        # finally the external functions if any
        for external_function in self._external_functions:

            path = external_function["function_path"]

            if external_function["is_composite"]:

                # we need to loop through the sub functions
                # can link them
                
                for i, primary_func in enumerate(new_model[path]._functions):

                    this_ef = external_function["external_functions"][i]

                    # for each function 
                    
                    if this_ef:

                        # if there are extrenal linked functions
                        
                        for fname, linked_function in this_ef.items():

                            # relink them
                            
                            primary_func.link_external_function(function=new_model[linked_function],
                                                                 internal_name=fname)
                        
            else:

                # do the same if it is not composite

                for fname, linked_function in external_function["external_functions"].items():

                    new_model[path].link_external_function(function=new_model[linked_function],
                                                            internal_name=fname)

                

                            

                    
        return new_model


class IndependentVariableParser(object):
    def __init__(self, name, definition):

        self._variable = parameter.IndependentVariable(name, **definition)

    def get_variable(self):

        return self._variable


class ParameterParser(object):
    def __init__(self, name, definition):

        self._links = []

        # NOTE: this is triggered only for parameters outside of functions

        if "prior" in definition:

            # Need the create a function for the prior first

            try:

                function_name = list(definition["prior"].keys())[0]
                parameters_definition = definition["prior"][function_name]

            except KeyError:  # pragma: no cover

                log.error("The prior for parameter %s is malformed" % name)

                raise ModelSyntaxError()

            # parse the function

            shape_parser = ShapeParser(name)

            prior_instance = shape_parser.parse(name, function_name,
                                                parameters_definition)

            # Substitute the definition with the instance, so that the following constructor will work
            definition["prior"] = prior_instance

        # Check if this is a linked parameter, i.e., if 'value' is something like f(source.spectrum.powerlaw.index)

        matches = re.findall("""f\((.+)\)""", str(definition["value"]))

        if matches:

            # This is an expression which marks a parameter
            # with a link to another parameter (or an IndependentVariable such as time)

            # Get the variable
            linked_variable = matches[0]

            # Now get the law

            if "law" not in definition:  # pragma: no cover

                log.error("The parameter %s in function %s "
                          " is linked to %s but lacks a 'law' attribute" %
                          (name, function_name, linked_variable))

                raise ModelSyntaxError()

            link_function_name = list(definition["law"].keys())[0]

            # ok, now we parse the linked parameter

            function_parser = ShapeParser(name)

            link_function_instance = function_parser.parse(
                name, link_function_name,
                definition["law"][link_function_name])

            self._links.append({
                "parameter_path": name,
                "law": link_function_instance,
                "variable": linked_variable,
            })

            # get rid of the 'law' entry

            definition.pop("law", None)

            # this parameter's value will be replaced later.
            # for now we just need to get rid of the f(param) entry

            definition["value"] = 1.0

        self._variable = parameter.Parameter(name, **definition)

    def get_variable(self):

        return self._variable

    @property
    def links(self):

        return self._links


class SourceParser(object):
    def __init__(self, source_name, source_definition):

        # Get the type of the source

        try:

            # Point source or extended source?

            source_type = re.findall(
                "\((%s|%s|%s)\)" % (SourceType.POINT_SOURCE, SourceType.EXTENDED_SOURCE, SourceType.PARTICLE_SOURCE),
                source_name,
            )[-1]

            
            
        except IndexError:  # pragma: no cover

            log.error("Don't recognize type for source '%s'. "
                "Valid types are '%s', '%s' or '%s'."
                % (source_name, SourceType.POINT_SOURCE, SourceType.EXTENDED_SOURCE, SourceType.PARTICLE_SOURCE))
            
            raise ModelSyntaxError()

        else:

            # Strip the source_type from the name

            source_name = source_name.split()[0]

        self._source_name = source_name

        # This will store the links (if any)
        self._links = []

        # This will store extra_setups (if any), used sometimes. For example, the function which uses naima
        # to make a synchrotron spectrum uses this to save and set up the particle distribution
        self._extra_setups = []

        # this will store any externally linked functions
        self._external_functions = []

        if source_type == SourceType.POINT_SOURCE.value:

            self._parsed_source = self._parse_point_source(source_definition)

        elif source_type == SourceType.EXTENDED_SOURCE.value:

            self._parsed_source = self._parse_extended_source(
                source_definition)

        elif source_type == SourceType.PARTICLE_SOURCE.value:

            self._parsed_source = self._parse_particle_source(
                source_definition)

    @property
    def extra_setups(self):

        return self._extra_setups

    @property
    def external_functions(self) -> List[Dict[str,str]]:

        return self._external_functions
    
    @property
    def links(self):

        return self._links

    def get_source(self):

        return self._parsed_source

    def _parse_particle_source(self, particle_source_definition):

        # Parse the spectral information

        try:

            spectrum = particle_source_definition["spectrum"]

        except KeyError:  # pragma: no cover

            log.error("Point source %s is missing the 'spectrum' attribute" %
                      self._source_name)

            raise ModelSyntaxError()

        components = []

        for component_name, component_definition in list(
                particle_source_definition["spectrum"].items()):

            this_component = self._parse_spectral_component(
                component_name, component_definition)

            components.append(this_component)

        this_particle_source = particle_source.ParticleSource(
            self._source_name, components=components)

        return this_particle_source

    def _parse_point_source(self, pts_source_definition):

        # Parse the positional information

        try:

            position_definition = pts_source_definition["position"]

        except KeyError:  # pragma: no cover

            log.error("Point source %s is missing the 'position' attribute" %
                      self._source_name)

            raise ModelSyntaxError()

        this_sky_direction = self._parse_sky_direction(position_definition)

        # Parse the spectral information

        try:

            spectrum = pts_source_definition["spectrum"]

        except KeyError:  # pragma: no cover

            log.error("Point source %s is missing the 'spectrum' attribute" %
                      self._source_name)

            raise ModelSyntaxError()

        components = []

        for component_name, component_definition in list(
                pts_source_definition["spectrum"].items()):

            try:

                this_component = self._parse_spectral_component(
                    component_name, component_definition)

                components.append(this_component)

            except:

                raise

        try:

            this_point_source = point_source.PointSource(
                self._source_name,
                sky_position=this_sky_direction,
                components=components,
            )

        except:

            raise

        return this_point_source

    def _parse_sky_direction(self, sky_direction_definition):

        # Instance the SkyDirection class using the coordinates provided

        coordinates = {}

        if "ra" in sky_direction_definition and "dec" in sky_direction_definition:

            par_parser = ParameterParser("ra", sky_direction_definition["ra"])

            ra = par_parser.get_variable()

            if ra.bounds == (None, None):
                ra.bounds = (0, 360)

            par_parser = ParameterParser("dec",
                                         sky_direction_definition["dec"])

            dec = par_parser.get_variable()

            if dec.bounds == (None, None):
                dec.bounds = (-90, 90)

            coordinates["ra"] = ra
            coordinates["dec"] = dec

        elif "l" in sky_direction_definition and "b" in sky_direction_definition:

            par_parser = ParameterParser("l", sky_direction_definition["l"])

            l = par_parser.get_variable()

            if l.bounds == (None, None):
                l.bounds = (0, 360)

            par_parser = ParameterParser("b", sky_direction_definition["b"])

            b = par_parser.get_variable()

            if b.bounds == (None, None):
                b.bounds = (-90, 90)

            coordinates["l"] = l
            coordinates["b"] = b

        else:  # pragma: no cover

            log.error("Position specification for source %s has an invalid coordinate pair. "
                " You need to specify either 'ra' and 'dec', or 'l' and 'b'."
                % self._source_name)
            
            raise ModelSyntaxError( )

        # Check if there is a equinox specification

        if "equinox" in sky_direction_definition:
            coordinates["equinox"] = sky_direction_definition["equinox"]

        try:

            this_sky_direction = sky_direction.SkyDirection(**coordinates)

        except sky_direction.WrongCoordinatePair:  # pragma: no cover

            log.error("Position specification for source %s has an invalid coordinate pair"
                % self._source_name)
            
            raise ModelSyntaxError( )

        return this_sky_direction

    def _parse_polarization(self, polarization_definititon):

        polarization_params = {}

        if "degree" in polarization_definititon and "angle" in polarization_definititon:

            par_parser = ParameterParser("degree", polarization_definititon["degree"])

            degree = par_parser.get_variable()

            degree.bounds = (0, 100)

            par_parser = ParameterParser("angle", polarization_definititon["angle"])

            angle = par_parser.get_variable()

            angle.bounds = (0, 180)

            this_polarization = polarization.LinearPolarization(
                angle=angle, degree=degree
            )

        elif (
            "I" in polarization_definititon
            and "U" in polarization_definititon
            and "Q" in polarization_definititon
            and "V" in polarization_definititon
        ):

            par_parser = ParameterParser("I", polarization_definititon["I"])

            I = par_parser.get_variable()

            I.bounds = (0, 1)

            par_parser = ParameterParser("U", polarization_definititon["U"])

            U = par_parser.get_variable()

            U.bounds = (0, 1)

            par_parser = ParameterParser("Q", polarization_definititon["Q"])

            Q = par_parser.get_variable()

            Q.bounds = (0, 1)

            par_parser = ParameterParser("V", polarization_definititon["V"])

            V = par_parser.get_variable()

            V.bounds = (0, 1)

            this_polarization = polarization.StokesPolarization(I=I, Q=Q, U=U, V=V)

        else:

            # just make a default polarization

            this_polarization = polarization.Polarization()
            # raise ModelSyntaxError("Polarization specification for source %s has an invalid parameters. "
            #                        " You need to specify either 'angle' and 'degree', or 'I' ,'Q', 'U' and 'V'."
            #                        % self._source_name)

        return this_polarization

    def _parse_spectral_component(self, component_name, component_definition):

        # Parse the shape definition, which is the first to occur

        try:

            function_name = list(component_definition.keys())[0]
            parameters_definition = component_definition[function_name]

        except KeyError:  # pragma: no cover

            log.error("The component %s of source %s is malformed"
                % (component_name, self._source_name))
            
            raise ModelSyntaxError( )

        # parse the function

        # now split the parameters and the properties

        
        shape_parser = ShapeParser(self._source_name)

        shape = shape_parser.parse(
            component_name, function_name, parameters_definition, is_spatial=False
        )

        # Get the links and extra setups, if any

        self._links.extend(shape_parser.links)
        self._extra_setups.extend(shape_parser.extra_setups)
        self._external_functions.extend(shape_parser.external_functions)

        
        if "polarization" in component_definition:

            # get the polarization

            polarization_definition = component_definition["polarization"]

            this_polarization = self._parse_polarization(polarization_definition)

        else:

            this_polarization = polarization.Polarization()

        this_spectral_component = spectral_component.SpectralComponent(
            component_name, shape, this_polarization
        )

        return this_spectral_component

    def _parse_extended_source(self, ext_source_definition):

        # The first item in the dictionary is the definition of the extended shape
        name_of_spatial_shape = list(ext_source_definition.keys())[0]

        spatial_shape_parser = ShapeParser(self._source_name)

        spatial_shape = spatial_shape_parser.parse(
            "n.a.",
            name_of_spatial_shape,
            list(ext_source_definition.values())[0],
            is_spatial=True,
        )

        # Get the links and extra setups, if any

        self._links.extend(spatial_shape_parser.links)
        self._extra_setups.extend(spatial_shape_parser.extra_setups)
        self._external_functions.extend(spatial_shape_parser.external_functions)
        # Parse the spectral information

        try:

            spectrum = ext_source_definition["spectrum"]

        except KeyError:  # pragma: no cover

            log.error("Ext. source %s is missing the 'spectrum' attribute" % self._source_name)
            
            raise ModelSyntaxError( )

        components = []

        for component_name, component_definition in list(
            ext_source_definition["spectrum"].items()
        ):
            this_component = self._parse_spectral_component(
                component_name, component_definition
            )

            components.append(this_component)

        this_ext_source = extended_source.ExtendedSource(
            self._source_name, spatial_shape, components=components
        )

        return this_ext_source


class ShapeParser(object):
    def __init__(self, source_name):

        self._source_name = source_name
        self._links = []
        self._extra_setups = []
        self._external_functions = []

        
    @property
    def links(self):

        return self._links

    @property
    def extra_setups(self):

        return self._extra_setups

    @property
    def external_functions(self):

        return self._external_functions
    
    def parse(
        self, component_name, function_name, parameters_definition, is_spatial=False
    ):

        return self._parse_shape_definition(
            component_name, function_name, parameters_definition, is_spatial
        )

    @staticmethod
    def _fix(value):
        # Remove new lines where it shouldn't be any
        # Sometimes YAML add new lines in the middle of definitions,
        # such as in units
        return value.replace("\n", " ")

    def _parse_shape_definition(
        self, component_name, function_name, parameters_definition, is_spatial=False
    ):

        # Get the function

        if "expression" in parameters_definition:

            # This is a composite function
            function_instance = function.get_function(
                function_name, parameters_definition["expression"]
            )

            is_composite = True
            
        else:

            try:

                function_instance = function.get_function(function_name)

                is_composite = False
                
            except function.UnknownFunction:  # pragma: no cover

                log.error( 
                    "Function %s, specified as shape for %s of source %s, is not a "
                    "known function"
                    % (function_name, component_name, self._source_name)
                )
                
                raise ModelSyntaxError()

        # Loop over the parameters of the function instance, instead of the specification,
        # so we can understand if there are parameters missing from the specification

        for parameter_name, _ in function_instance.parameters.items():

            try:

                this_definition = parameters_definition[parameter_name]

            except KeyError:  # pragma: no cover

                log.error( 
                    "Function %s, specified as shape for %s of source %s, lacks "
                    "the definition for parameter %s"
                    % (function_name, component_name, self._source_name, parameter_name)
                )

                for k,v in parameters_definition.items():
                
                    log.error((k, v))
                
                raise ModelSyntaxError()

            # Update the parameter. Note that the order is important, because trying to set the value before the
            # minimum and maximum could result in a error.

            # All these specifications are optional. If they are not present, then the default value
            # already contained in the instance of the function will be used

            # Ignore for a second the RuntimeWarning that is printed if the default value in the function definition
            # is outside the bounds defined here

            with warnings.catch_warnings():

                warnings.simplefilter("ignore", RuntimeWarning)

                if "min_value" in this_definition:
                    function_instance.parameters[
                        parameter_name
                    ].min_value = this_definition["min_value"]

                if "max_value" in this_definition:
                    function_instance.parameters[
                        parameter_name
                    ].max_value = this_definition["max_value"]

            if "delta" in this_definition:
                function_instance.parameters[parameter_name].delta = this_definition[
                    "delta"
                ]

            if "free" in this_definition:
                function_instance.parameters[parameter_name].free = this_definition[
                    "free"
                ]

            if "unit" in this_definition:
                function_instance.parameters[parameter_name].unit = self._fix(
                    this_definition["unit"]
                )

            # Now set the value, which must be present

            if "value" not in this_definition:  # pragma: no cover


                log.error( 
                    "The parameter %s in function %s, specified as shape for %s "
                    "of source %s, lacks a 'value' attribute"
                    % (parameter_name, function_name, component_name, self._source_name))
                
                raise ModelSyntaxError( 
                )

            # Check if this is a linked parameter, i.e., if 'value' is something like f(source.spectrum.powerlaw.index)

            matches = re.findall("""f\((.+)\)""", str(this_definition["value"]))

            if matches:

                # This is an expression which marks a parameter
                # with a link to another parameter (or an IndependentVariable such as time)

                # Get the variable
                linked_variable = matches[0]

                # Now get the law

                if "law" not in this_definition:  # pragma: no cover

                    log.error( 
                        "The parameter %s in function %s, specified as shape for %s "
                        "of source %s, is linked to %s but lacks a 'law' attribute"
                        % (
                            parameter_name,
                            function_name,
                            component_name,
                            self._source_name,
                            linked_variable,
                        )
                    )
                    
                    raise ModelSyntaxError()

                link_function_name = list(this_definition["law"].keys())[0]

                link_function_instance = self._parse_shape_definition(
                    component_name,
                    link_function_name,
                    this_definition["law"][link_function_name],
                )

                if is_spatial:
                    path = ".".join([self._source_name, function_name, parameter_name])
                else:
                    path = ".".join(
                        [
                            self._source_name,
                            "spectrum",
                            component_name,
                            function_name,
                            parameter_name,
                        ]
                    )

                self._links.append(
                    {
                        "parameter_path": path,
                        "law": link_function_instance,
                        "variable": linked_variable,
                    }
                )

            else:

                # This is a normal (not linked) parameter

                function_instance.parameters[parameter_name].value = this_definition[
                    "value"
                ]

            # Setup the prior for this parameter, if it exists
            if "prior" in this_definition:
                # Get the function for this prior

                # A name to display in case of errors

                name_for_errors = (
                    "prior for %s" % function_instance.parameters[parameter_name].path
                )

                prior_function_name = list(this_definition["prior"].keys())[0]

                prior_function_definition = this_definition["prior"][
                    prior_function_name
                ]

                prior_function = self._parse_shape_definition(
                    name_for_errors, prior_function_name, prior_function_definition
                )

                # Set it as prior for current parameter

                function_instance.parameters[parameter_name].prior = prior_function

        
        
                
        if function_instance.has_properties:

            # now collect the properties
            
            # the properties are stored in the parameters defintion
            # as well
            
            for property_name, _ in function_instance.properties.items():

                try:

                    this_definition = parameters_definition[property_name]

                except KeyError:  # pragma: no cover

                    log.error( 
                        "Function %s, specified as shape for %s of source %s, lacks "
                        "the definition for property %s"
                        % (function_name, component_name, self._source_name, property_name)
                    )

                    for k,v in parameters_definition.items():

                        log.error((k, v))

                    raise ModelSyntaxError()


                if "value" not in this_definition:

                    log.error("The property %s in function %s, specified as shape for %s "
                              "of source %s, lacks a 'value' attribute"
                    % (property_name, function_name, component_name, self._source_name))
                
                    raise ModelSyntaxError()


                function_instance.properties[property_name].value = this_definition['value']
                    
                

                
        # Now handle extra_setup if any
        if "extra_setup" in parameters_definition:

            if is_spatial:
                path = ".".join([self._source_name, function_name])
            else:
                path = ".".join(
                    [self._source_name, "spectrum", component_name, function_name]
                )

            self._extra_setups.append(
                {
                    "function_path": path,
                    "extra_setup": parameters_definition["extra_setup"],
                }
            )
        if "external_functions" in parameters_definition:

            if is_spatial:
                path = ".".join([self._source_name, function_name])
            else:
                path = ".".join(
                    [self._source_name, "spectrum", component_name, function_name])


            
            self._external_functions.append(
                {"function_path": path,                 
                 "external_functions": parameters_definition["external_functions"],
                 "is_composite": is_composite
            }
            )

            
        return function_instance
