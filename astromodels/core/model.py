__author__ = 'giacomov'

import collections

import os
import pandas as pd
import warnings

from astromodels.core.my_yaml import my_yaml
from astromodels.core.parameter import Parameter, IndependentVariable
from astromodels.core.tree import Node, DuplicatedNode
from astromodels.functions.function import get_function
from astromodels.sources.source import Source, POINT_SOURCE, EXTENDED_SOURCE, PARTICLE_SOURCE
from astromodels.utils.disk_usage import disk_usage
from astromodels.utils.long_path_formatter import long_path_formatter


class ModelFileExists(IOError):
    pass


class InvalidInput(ValueError):
    pass


class CannotWriteModel(IOError):
    def __init__(self, directory, message):
        # Add a report on disk usage to the message

        free_space = disk_usage(directory).free

        message += "\nFree space on the file system hosting %s was %.2f Mbytes" % (
            directory, free_space / 1024.0 / 1024.0)

        super(CannotWriteModel, self).__init__(message)


class ModelInternalError(ValueError):

    pass


class Model(Node):

    def __init__(self, *sources):

        # There must be at least one source
        assert len(sources) > 0, "You need to have at least one source in the model."

        # Setup the node, using the special name '__root__' to indicate that this is the root of the tree

        super(Model, self).__init__("__root__")

        # Dictionary to keep point sources

        self._point_sources = collections.OrderedDict()

        # Dictionary to keep extended sources

        self._extended_sources = collections.OrderedDict()

        # Dictionary to keep particle sources

        self._particle_sources = collections.OrderedDict()

        # Loop over the provided sources and process them

        for source in sources:

            self._add_source(source)

        # Now make the list of all the existing parameters

        self._update_parameters()

        # This controls the verbosity of the display
        self._complete_display = False

    def _add_source(self, source):
        """
        Remember to call _update_parameters after this!
        :param source:
        :return:
        """

        try:

            self._add_child(source)

        except AttributeError:

            if isinstance(source, Source):

                raise DuplicatedNode("More than one source with the name '%s'. You cannot use the same name for multiple "
                                     "sources" % source.name)

            else:  # pragma: no cover

                raise

        # Now see if this is a point or extended source, and add them to the
        # appropriate dictionary

        if source.source_type == POINT_SOURCE:

            self._point_sources[source.name] = source

        elif source.source_type == EXTENDED_SOURCE:

            self._extended_sources[source.name] = source

        elif source.source_type == PARTICLE_SOURCE:

            self._particle_sources[source.name] = source

        else:  # pragma: no cover

            raise InvalidInput("Input sources must be either a point source or an extended source")

    def _remove_source(self, source_name):
        """
        Remember to call _update_parameters after this
        :param source_name:
        :return:
        """

        assert source_name in self.sources, "Source %s is not part of the current model" % source_name

        source = self.sources.pop(source_name)

        if source.source_type == POINT_SOURCE:

            self._point_sources.pop(source.name)

        elif source.source_type == EXTENDED_SOURCE:

            self._extended_sources.pop(source.name)

        elif source.source_type == PARTICLE_SOURCE:

            self._particle_sources.pop(source.name)

        self._remove_child(source_name)

    def _find_parameters(self, node):

        instances = collections.OrderedDict()

        for child in node._get_children():

            if isinstance(child, Parameter):

                path = child._get_path()

                instances[path] = child

                for sub_child in child._get_children():

                    instances.update(self._find_parameters(sub_child))

            else:

                instances.update(self._find_parameters(child))

        return instances

    def _update_parameters(self):

        self._parameters = self._find_parameters(self)

    @property
    def parameters(self):
        """
        Return a dictionary with all parameters

        :return: dictionary of parameters
        """
        self._update_parameters()

        return self._parameters

    @property
    def free_parameters(self):
        """
        Get a dictionary with all the free parameters in this model

        :return: dictionary of free parameters
        """

        # Refresh the list

        self._update_parameters()

        # Filter selecting only free parameters

        free_parameters_dictionary = collections.OrderedDict()

        for parameter_name, parameter in self._parameters.iteritems():

            if parameter.free:

                free_parameters_dictionary[parameter_name] = parameter

        return free_parameters_dictionary

    @property
    def linked_parameters(self):
        """
        Get a dictionary with all parameters in this model in a linked status. A parameter is in a linked status
        if it is linked to another parameter (i.e. it is forced to have the same value of the other parameter), or
        if it is linked with another parameter or an independent variable through a law.

        :return: dictionary of linked parameters
        """

        # Refresh the list

        self._update_parameters()

        # Filter selecting only free parameters

        linked_parameter_dictionary = collections.OrderedDict()

        for parameter_name, parameter in self._parameters.iteritems():

            if parameter.has_auxiliary_variable():

                linked_parameter_dictionary[parameter_name] = parameter

        return linked_parameter_dictionary

    def set_free_parameters(self, values):
        """
        Set the free parameters in the model to the provided values.

        NOTE: of course, order matters

        :param values: a list of new values
        :return: None
        """

        assert len(values) == len(self.free_parameters)

        for parameter, this_value in zip(self.free_parameters.values(), values):

            parameter.value = this_value

    def __getitem__(self, path):
        """
        Get a parameter from a path like "source_1.component.powerlaw.logK". This might be useful in certain
        context, although in an interactive analysis there is no reason to use this.

        :param path: the address of the parameter
        :return: the parameter
        """

        return self._get_child_from_path(path)

    def __contains__(self, path):
        """
        This allows the model to be used with the "in" operator, like;

        > if 'myparameter' in model:
        >    print("Myparameter is contained in the model")

        :param path: the parameter to look for
        :return:
        """

        try:

            _ = self._get_child_from_path(path)

        except (AttributeError, KeyError):

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
    def point_sources(self):
        """
        Returns the dictionary of all defined point sources

        :return: collections.OrderedDict()
        """
        return self._point_sources

    @property
    def extended_sources(self):
        """
        Returns the dictionary of all defined extended sources

        :return: collections.OrderedDict()

        """
        return self._extended_sources

    @property
    def particle_sources(self):
        """
        Returns the dictionary of all defined particle sources

        :return: collections.OrderedDict()

        """
        return self._particle_sources

    @property
    def sources(self):
        """
        Returns a dictionary containing all defined sources (of any kind)

        :return: collections.OrderedDict()

        """

        sources = collections.OrderedDict()

        for d in (self.point_sources, self.extended_sources, self.particle_sources):

            sources.update(d)

        return sources

    def add_source(self, new_source):
        """
        Add the provided source to the model

        :param new_source: the new source to be added (an instance of PointSource, ExtendedSource or ParticleSource)
        :return: (none)
        """

        self._add_source(new_source)

        self._update_parameters()

    def remove_source(self, source_name):
        """
        Returns a new model with the provided source removed from the current model

        :param source_name: the name of the source to be removed
        :return: a new Model instance without the source
        """

        self._remove_source(source_name)

        self._update_parameters()

    def add_independent_variable(self, variable):
        """
        Add a global independent variable to this model, such as time.

        :param variable: an IndependentVariable instance
        :return: none
        """

        assert isinstance(variable, IndependentVariable),"Variable must be an instance of IndependentVariable"

        if self._has_child(variable.name):

            self._remove_child(variable.name)

        self._add_child(variable)

    def remove_independent_variable(self, variable_name):
        """
        Remove an independent variable which was added with add_independent_variable

        :param variable_name: name of variable to remove
        :return:
        """

        self._remove_child(variable_name)

    def add_external_parameter(self, parameter):
        """
        Add a parameter that comes from something other than a function, to the model.

        :param parameter: a Parameter instance
        :return: none
        """

        assert isinstance(parameter, Parameter), "Variable must be an instance of IndependentVariable"

        if self._has_child(parameter.name):

            # Remove it from the children only if it is a Parameter instance, otherwise don't, which will
            # make the _add_child call fail (which is the expected behaviour! You shouldn't call two children
            # with the same name)

            if isinstance(self._get_child(parameter.name), Parameter):

                warnings.warn("External parameter %s already exist in the model. Overwriting it..." % parameter.name,
                              RuntimeWarning)

                self._remove_child(parameter.name)

        # This will fail if another node with the same name is already in the model

        self._add_child(parameter)

    def remove_external_parameter(self, parameter_name):
        """
        Remove an external parameter which was added with add_external_parameter

        :param variable_name: name of parameter to remove
        :return:
        """

        self._remove_child(parameter_name)

    def link(self, parameter_1, parameter_2, link_function=None):
        """
        Link the value of the provided parameters through the provided function (identity is the default, i.e.,
        parameter_1 = parameter_2).

        :param parameter_1: the first parameter
        :param parameter_2: the second parameter
        :param link_function: a function instance. If not provided, the identity function will be used by default.
        Otherwise, this link will be set: parameter_1 = link_function(parameter_2)
        :return: (none)
        """

        if link_function is None:
            # Use the Line function by default, with both parameters fixed so that the two
            # parameters to be linked will vary together
            link_function = get_function('Line')

            link_function.a.value = 1
            link_function.a.fix = True

            link_function.b.value = 0
            link_function.b.fix = True

        parameter_1.add_auxiliary_variable(parameter_2, link_function)

        # Now set the units of the link function
        link_function.set_units(parameter_2.unit, parameter_1.unit)

    def unlink(self, parameter):
        """
        Sets free a parameter which has been linked previously

        :param parameter: the parameter to be set free
        :return: (none)
        """

        if parameter.has_auxiliary_variable():

            parameter.remove_auxiliary_variable()

        else:

            with warnings.catch_warnings():

                warnings.simplefilter("always", RuntimeWarning)

                warnings.warn("Parameter %s has no link to be removed." % parameter.path, RuntimeWarning)

    def display(self, complete=False):
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

            new_line = '<br>'

        else:

            new_line = '\n'

        # Table with the summary of the various kind of sources
        sources_summary = pd.DataFrame.from_items((('Point sources', [self.get_number_of_point_sources()]),
                                                   ('Extended sources', [self.get_number_of_extended_sources()]),
                                                   ('Particle sources', [self.get_number_of_particle_sources()])),
                                                  columns=['N'], orient='index')

        # These properties traverse the whole tree everytime, so let's cache their results here
        parameters = self.parameters
        free_parameters = self.free_parameters
        linked_parameters = self.linked_parameters

        # Summary of free parameters
        if len(free_parameters) > 0:

            parameter_dict = collections.OrderedDict()

            for parameter_name, parameter in free_parameters.iteritems():
                # Generate table with only a minimal set of info

                # Generate table with only a minimal set of info
                if rich_output:

                    this_name = long_path_formatter(parameter_name, 70)

                else:

                    # In a terminal we need to use less characters

                    this_name = long_path_formatter(parameter_name, 40)

                d = parameter.to_dict()
                parameter_dict[this_name] = collections.OrderedDict()

                for key in ['value','unit','min_value','max_value']:

                    parameter_dict[this_name][key] = d[key]

            free_parameters_summary = pd.DataFrame.from_dict(parameter_dict).T

            # Re-order it
            free_parameters_summary = free_parameters_summary[['value','min_value','max_value','unit']]

        else:

            free_parameters_summary = pd.DataFrame()

        if len(parameters) - len(free_parameters) - len(linked_parameters) > 0:

            fixed_parameter_dict = collections.OrderedDict()

            for parameter_name, parameter in parameters.iteritems():

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

                for key in ['value','unit','min_value','max_value']:

                    fixed_parameter_dict[this_name][key] = d[key]

            fixed_parameters_summary = pd.DataFrame.from_dict(fixed_parameter_dict).T

            # Re-order it
            fixed_parameters_summary = fixed_parameters_summary[['value','min_value','max_value','unit']]

        else:

            fixed_parameters_summary = pd.DataFrame()

        # Summary of linked parameters

        linked_frames = []

        if linked_parameters:

            for parameter_name, parameter in linked_parameters.iteritems():

                parameter_dict = collections.OrderedDict()

                # Generate table with only a minimal set of info

                variable, law = parameter.auxiliary_variable

                this_dict = collections.OrderedDict()

                this_dict['linked to'] = variable.path
                this_dict['function'] = law.name
                this_dict['current value'] = parameter.value
                this_dict['unit'] = parameter.unit

                parameter_dict[parameter_name] = this_dict

                this_parameter_frame = pd.DataFrame.from_dict(parameter_dict)

                linked_frames.append(this_parameter_frame)

        else:

            # No linked parameters

            pass

        empty_frame = "(none)%s" % new_line

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

            if len(linked_frames) == 0:

                linked_summary_representation = empty_frame

            else:

                linked_summary_representation = ""

                for linked_frame in linked_frames:

                    linked_summary_representation += linked_frame.__repr__()
                    linked_summary_representation += "%s%s" % (new_line, new_line)

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

        representation += "%sFree parameters (%i):%s" % (new_line, len(free_parameters), new_line)

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

            representation += '(abridged. Use complete=True to see all fixed parameters)%s' % new_line

        representation += new_line

        # Linked parameters

        representation += "%sLinked parameters (%i):%s" % (new_line, len(self.linked_parameters), new_line)

        if not rich_output:

            representation += "----------------------%s%s" % (new_line, new_line)

        else:

            representation += new_line

        representation += linked_summary_representation

        return representation

    def to_dict_with_types(self):

        # Get the serialization dictionary

        data = self.to_dict()

        # Add the types to the sources

        for key in data.keys():

            try:

                element = self._get_child(key)

            except KeyError:  # pragma: no cover

                raise RuntimeError("Source %s is unknown" % key)

            else:

                # There are three possible cases. Either the element is a source, or it is an independent
                # variable, or a parameter

                if hasattr(element, 'source_type'):

                    # Change the name of the key adding the source type

                    data['%s (%s)' % (key, element.source_type)] = data.pop(key)

                elif isinstance(element, IndependentVariable):

                    data['%s (%s)' % (key, 'IndependentVariable')] = data.pop(key)

                elif isinstance(element, Parameter):

                    data['%s (%s)' % (key, 'Parameter')] = data.pop(key)

                else: # pragma: no cover

                    raise ModelInternalError("Found an unknown class at the top level")

        return data

    def save(self, output_file, overwrite=False):

        """Save the model to disk"""

        if os.path.exists(output_file) and overwrite is False:

            raise ModelFileExists("The file %s exists already. If you want to overwrite it, use the 'overwrite=True' "
                                  "options as 'model.save(\"%s\", overwrite=True)'. " % (output_file, output_file))

        else:

            data = self.to_dict_with_types()

            # Write it to disk

            try:

                # Get the YAML representation of the data

                representation = my_yaml.dump(data)

                with open(output_file, "w+") as f:

                    # Add a new line at the end of each voice (just for clarity)

                    f.write(representation.replace("\n", "\n\n"))

            except IOError:

                raise CannotWriteModel(os.path.dirname(os.path.abspath(output_file)),
                                       "Could not write model file %s. Check your permissions to write or the "
                                       "report on the free space which follows: " % output_file)

    def get_number_of_point_sources(self):
        """
        Return the number of point sources

        :return: number of point sources
        """
        return len(self._point_sources)

    def get_point_source_position(self, id):
        """
        Get the point source position (R.A., Dec)

        :param id: id of the source
        :return: a tuple with R.A. and Dec.
        """

        pts = self._point_sources.values()[id]

        return pts.position.get_ra(), pts.position.get_dec()

    def get_point_source_fluxes(self, id, energies):
        """
        Get the fluxes from the id-th point source

        :param id: id of the source
        :param energies: energies at which you need the flux
        :return: fluxes
        """

        return self._point_sources.values()[id](energies)

    def get_point_source_name(self, id):

        return self._point_sources.values()[id].name

    def get_number_of_extended_sources(self):
        """
        Return the number of extended sources

        :return: number of extended sources
        """
        return len(self._extended_sources)

    def get_extended_source_fluxes(self, id, j2000_ra, j2000_dec, energies):
        """
        Get the flux of the id-th extended sources at the given position at the given energies

        :param id: id of the source
        :param j2000_ra: R.A. where the flux is desired
        :param j2000_dec: Dec. where the flux is desired
        :param energies: energies at which the flux is desired
        :return: flux array
        """

        return self._extended_sources.values()[id](j2000_ra, j2000_dec, energies)

    def get_extended_source_name(self, id):
        """
        Return the name of the n-th extended source

        :param id: id of the source (integer)
        :return: the name of the id-th source
        """

        return self._extended_sources.values()[id].name

    def get_extended_source_boundaries(self, id):

        (ra_min, ra_max), (dec_min, dec_max) = self._extended_sources.values()[id].get_boundaries()

        return ra_min, ra_max, dec_min, dec_max

    def is_inside_any_extended_source(self, j2000_ra, j2000_dec):

        for ext_source in self.extended_sources.values():

            (ra_min, ra_max), (dec_min, dec_max) = ext_source.get_boundaries()

            # Transform from the 0...360 convention to the -180..180 convention, so that
            # the comparison is easier
            if ra_min > 180:

                ra_min = -(360-ra_min)

            if ra_min <= j2000_ra <= ra_max and dec_min <= j2000_dec <= dec_max:

                return True

        # If we are here, it means that no extended source contains the provided coordinates

        return False

    def get_number_of_particle_sources(self):
        """
        Return the number of particle sources

        :return: number of particle sources
        """

        return len(self._particle_sources)

    def get_particle_source_fluxes(self, id, energies):
        """
        Get the fluxes from the id-th point source

        :param id: id of the source
        :param energies: energies at which you need the flux
        :return: fluxes
        """

        return self._particle_sources.values()[id](energies)

    def get_particle_source_name(self, id):

        return self._particle_sources.values()[id].name