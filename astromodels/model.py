__author__ = 'giacomov'

import collections
import os

from astromodels.dual_access_class import DualAccessClass
from astromodels.sources.point_source import PointSource
from astromodels.sources.extended_source import ExtendedSource
from astromodels.my_yaml import my_yaml
from astromodels.utils.disk_usage import disk_usage
from astromodels.utils.table import dict_to_table
from astromodels.utils.io import display


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


class Model(DualAccessClass):
    def __init__(self, *sources):

        # Dictionary to keep all parameters

        self._parameters = collections.OrderedDict()

        # Dictionary to keep point sources

        self._point_sources = collections.OrderedDict()

        # Dictionary to keep extended sources

        self._extended_sources = collections.OrderedDict()

        # This is needed to allow the user to access the sources by name as model.source, independently of their
        # type

        self._sources = collections.OrderedDict()

        # Loop over the provided sources and process them

        for source in sources:

            assert source.name not in self._sources, "You cannot use the same name for multiple sources."

            self._sources[source.name] = source

            # Now see if this is a point or extended source, and add them to the
            # appropriate dictionary

            if isinstance(source, PointSource):

                self._point_sources[source.name] = source

            elif isinstance(source, ExtendedSource):

                self._extended_sources[source.name] = source

            else:

                raise InvalidInput("Input sources must be instances of either PointSource or ExtendedSource")

        # Also create the lists for the point and extended sources

        self._point_sources_list = self._point_sources.values()
        self._extended_sources_list = self._extended_sources.values()

        super(Model, self).__init__('source', self._sources)

    @property
    def free_parameters(self):
        """
        Get a dictionary with all the free parameters in this model

        :return: dictionary of free parameters
        """

        parameters_dictionary = self._user_friendly_repr_of_parameters(self.to_dict())

        # Filter out fixed parameters

        free_parameters_dictionary = collections.OrderedDict()

        for param_name, param_dict in parameters_dictionary.iteritems():

            if param_dict['free']:
                param = self.get_parameter(param_name)

                free_parameters_dictionary[param_name] = param

        return free_parameters_dictionary

    def get_parameter(self, user_friendly_parameter_address):
        """
        Get a parameter from an address like "source_1.component.powerlaw.logK". This might be useful in certain
        context, although in an interactive analysis there is no reason to use this.

        :param user_friendly_parameter_address: the address of the parameter
        :return: the parameter
        """

        return eval('self.%s' % user_friendly_parameter_address)

    @property
    def point_sources(self):
        return self._point_sources

    @property
    def extended_sources(self):
        return self._extended_sources

    @property
    def sources(self):
        return self._sources

    def __repr__base(self, rich_output=False):

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

        representation += "%s%sParameters:%s" % (div, div, div)

        # Divide the dictionary source by source, to make the display clearer

        dict_representation = self.to_dict()

        for source_name, source in dict_representation.iteritems():

            representation += "%s%s%s" % (div, source_name, div)

            # Create a dictionary with only one key

            this_dict = collections.OrderedDict()
            this_dict[source_name] = source

            param_dict = self._user_friendly_repr_of_parameters(this_dict, add_source_name=False)

            table = dict_to_table(param_dict)

            if rich_output:

                representation += table._repr_html_()

            else:

                representation += table.__repr__()

            representation += div

        return representation

    def __repr__(self):

        return self.__repr__base(rich_output=False)

    def _repr_html_(self):

        return self.__repr__base(rich_output=True)

    def display(self):

        display(self)

    def to_dict(self):
        """Prepare dictionary for serialization"""

        # Prepare data

        data = collections.OrderedDict()

        # Loop over all sources and get the serializations to dictionaries

        for source_name, source in self.sources.iteritems():
            data[source_name] = source.to_dict()

        return data

    def save(self, output_file, overwrite=False):

        """Save the model to disk"""

        if os.path.exists(output_file) and overwrite is False:

            raise ModelFileExists("The file %s exists already. If you want to overwrite it, use the 'overwrite=True' "
                                  "options as 'model.save(\"%s\", overwrite=True)'. " % (output_file, output_file))

        else:

            # Get the serialization dictionary

            data = self.to_dict()

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

    @staticmethod
    def _user_friendly_repr_of_parameters(sources, add_source_name=True):
        """
        Returns a dictionary with a user-friendly key/value representation of the model.

        :return: dictionary of parameters with keys which can be directly used by the user to access parameters
        """

        # NOTE: there are things in the definition of a source obtained with to_dict() which might
        # not be Parameter instances. For example, the equinox in a sky direction or the file name
        # in an extended source.

        params = collections.OrderedDict()

        for src_name, src in sources.iteritems():

            # Get the type of the source

            src_type = src.keys()[0]

            # Add the source name to the final key, unless add_source_name is False

            if add_source_name:

                master_key = "%s." % src_name

            else:

                master_key = ''

            # Check whether this is a point source or an extended one

            if src_type == 'point source':

                # Get the position parameters

                for param_name, param in src['point source']['position'].iteritems():

                    # Check whether this is effectively a parameter

                    if isinstance(param, dict):
                        # Generate the key and add it to the dictionary

                        key = '%sposition.%s' % (master_key, param_name)

                        params[key] = param

                # Loop over the components

                for component_name, component in src['point source']['spectrum'].iteritems():

                    # Spectral function

                    function_name = component['shape'].keys()[0]

                    function_parameters = component['shape'].values()[0]

                    for param_name, param in function_parameters.iteritems():

                        if isinstance(param, dict):
                            key = '%s%s.%s.%s' % (master_key, component_name, function_name, param_name)

                            params[key] = param

            else:

                raise NotImplemented("Extended sources not implemented")

        return params

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

        return pts.position.ra.value, pts.position.dec.value

    def get_point_source_fluxes(self, id, energies):
        """
        Get the fluxes from the id-th point source

        :param id: id of the source
        :param energies: energies at which you need the flux
        :return: fluxes
        """

        return self._point_sources_list[id].get_flux(energies)

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

        return self._extended_sources_list[id].get_brightness(j2000_ra, j2000_dec, energies)

    def get_extended_source_name(self, id):
        """
        Return the name of the n-th extended source

        :param id: id of the source (integer)
        :return: the name of the id-th source
        """

        return self._extended_sources_list[id]

    def get_extended_source_boundaries(self, id):

        ra_min, ra_max, dec_min, dec_max = self._extended_sources_list[id].get_boundaries()

        return ra_min, ra_max, dec_min, dec_max

    def is_inside_any_extended_source(self, j2000_ra, j2000_dec):

        return True
