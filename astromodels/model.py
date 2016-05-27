__author__ = 'giacomov'

import collections
import os
import warnings

from astromodels.sources.source import POINT_SOURCE, EXTENDED_SOURCE, PARTICLE_SOURCE

from astromodels.my_yaml import my_yaml
from astromodels.utils.disk_usage import disk_usage
from astromodels.utils.table import dict_to_table
from astromodels.parameter import Parameter, IndependentVariable
from astromodels.tree import Node
from astromodels.functions.function import get_function


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


class Model(Node):

    def __init__(self, *sources):

        # Setup the node

        super(Model, self).__init__("root")

        # Dictionary to keep point sources

        self._point_sources = collections.OrderedDict()

        # Dictionary to keep extended sources

        self._extended_sources = collections.OrderedDict()

        # Dictionary to keep particle sources

        self._particle_sources = collections.OrderedDict()

        # Loop over the provided sources and process them

        for source in sources:

            self._add_child(source)

            # Now see if this is a point or extended source, and add them to the
            # appropriate dictionary

            if source.source_type == POINT_SOURCE:

                self._point_sources[source.name] = source

            elif source.source_type == EXTENDED_SOURCE:

                self._extended_sources[source.name] = source

            elif source.source_type == PARTICLE_SOURCE:

                self._particle_sources[source.name] = source

            else:

                raise InvalidInput("Input sources must be either a point source or an extended source")

        # Also create the lists for the point and extended sources, for speed

        self._point_sources_list = self._point_sources.values()
        self._extended_sources_list = self._extended_sources.values()
        self._particle_sources_list = self._particle_sources.values()

        # Now make the list of all the existing parameters

        self._update_parameters()

    def _update_parameters(self):

        self._parameters = self._find_instances(Parameter)

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

    @property
    def point_sources(self):
        return self._point_sources

    @property
    def extended_sources(self):
        return self._extended_sources

    @property
    def particle_sources(self):
        return self._particle_sources

    def add_independent_variable(self, variable):
        """
        Add a global independent variable to this model, such as time.

        :param variable: an IndependentVariable instance
        :return: none
        """

        assert isinstance(variable, IndependentVariable),"Variable must be an instance of IndependentVariable"

        if variable.name in self._children:

            self._remove_child(variable.name)

        self._add_child(variable)

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
            # Use the identity function by default

            link_function = get_function('identity')

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

            warnings.warn("Parameter %s has no link to be removed." % parameter.path)

    def _repr__base(self, rich_output=False):

        if rich_output:

            div = '<br>'

        else:

            div = '\n'

        # Print the name of point sources if there are any

        representation = "Point sources: "

        if self._point_sources:

            representation += ",".join(self._point_sources.keys())

        else:

            representation += "(none)"

        # Print the name of extended sources if there are any

        representation += "%s%sExtended sources: " % (div, div)

        if self._extended_sources:

            representation += ",".join(self._extended_sources.keys())

        else:

            representation += "(none)"

        # Print the name of extended sources if there are any

        representation += "%s%sParticle sources: " % (div, div)

        if self._particle_sources:

            representation += ",".join(self._particle_sources.keys())

        else:

            representation += "(none)"

        representation += "%s%sFree parameters:%s" % (div, div, div)

        parameters_dict = collections.OrderedDict()

        for parameter_name, parameter in self.free_parameters.iteritems():

            # Generate table with only a minimal set of info

            parameters_dict[parameter_name] = parameter.to_dict()

        table = dict_to_table(parameters_dict, ['value','unit','min_value','max_value','delta','free'])

        if rich_output:

            representation += table._repr_html_()

        else:

            representation += table.__repr__()

            representation += div

        # Print linked parameters

        linked_parameters = self.linked_parameters

        if linked_parameters:

            representation += "%s%sLinked parameters:%s" % (div, div, div)

            parameters_dict = collections.OrderedDict()

            for parameter_name, parameter in linked_parameters.iteritems():

                # Generate table with only a minimal set of info

                variable, law = parameter.auxiliary_variable

                this_dict = collections.OrderedDict()

                this_dict['linked to'] = variable.path
                this_dict['function'] = law.name
                this_dict['current value'] = parameter.value.value
                this_dict['unit'] = parameter.unit

                parameters_dict[parameter_name] = this_dict

            table = dict_to_table(parameters_dict)

            if rich_output:

                representation += table._repr_html_()

            else:

                representation += table.__repr__()

                representation += div

        return representation

    def save(self, output_file, overwrite=False):

        """Save the model to disk"""

        if os.path.exists(output_file) and overwrite is False:

            raise ModelFileExists("The file %s exists already. If you want to overwrite it, use the 'overwrite=True' "
                                  "options as 'model.save(\"%s\", overwrite=True)'. " % (output_file, output_file))

        else:

            # Get the serialization dictionary

            data = self.to_dict()

            # Add the types to the sources

            for key in data.keys():

                try:

                    element = self._children[key]

                except KeyError:

                    raise RuntimeError("Source %s is unknown" % key)

                else:

                    # There are two possible cases. Either the element is a source, or it is an independent
                    # variable

                    if hasattr(element, 'source_type'):

                        # Change the name of the key adding the source type

                        data['%s (%s)' % (key, element.source_type)] = data.pop(key)

                    elif isinstance(element, IndependentVariable):

                        data['%s (%s)' % (key, 'IndependentVariable')] = data.pop(key)

                    else:

                        raise CannotWriteModel("Found an unknown class at the top level")

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

        pts = self._point_sources_list[id]

        return pts.position.get_ra(), pts.position.get_dec()

    def get_point_source_fluxes(self, id, energies):
        """
        Get the fluxes from the id-th point source

        :param id: id of the source
        :param energies: energies at which you need the flux
        :return: fluxes
        """

        return self._point_sources_list[id](energies)

    def get_point_source_name(self, id):

        return self._point_sources_list[id].name

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

        return self._extended_sources_list[id](j2000_ra, j2000_dec, energies)

    def get_extended_source_name(self, id):
        """
        Return the name of the n-th extended source

        :param id: id of the source (integer)
        :return: the name of the id-th source
        """

        return self._extended_sources_list[id]

    def get_extended_source_boundaries(self, id):

        (ra_min, ra_max), (dec_min, dec_max) = self._extended_sources_list[id].get_boundaries()

        return ra_min, ra_max, dec_min, dec_max

    def is_inside_any_extended_source(self, j2000_ra, j2000_dec):

        return True

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

        return self._particle_sources_list[id]._get_flux(energies)

    def get_particle_source_name(self, id):

        return self._particle_sources_list[id].name