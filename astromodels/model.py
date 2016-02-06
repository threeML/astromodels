__author__ = 'giacomov'

import collections
import os
from astromodels.dual_access_class import DualAccessClass
from astromodels.sources.point_source import PointSource
from astromodels.sources.extended_source import ExtendedSource
from astromodels.my_yaml import my_yaml
from astromodels.utils.disk_usage import disk_usage


class ModelFileExists(IOError):
    pass


class CannotWriteModel(IOError):

    def __init__(self, directory, message):

        # Add a report on disk usage to the message

        free_space = disk_usage(directory).free

        message += "\nFree space on the file system hosting %s was %.2f Mbytes" % (directory, free_space/1024.0/1024.0)

        super(CannotWriteModel, self).__init__(message)


class Model(DualAccessClass):

    def __init__(self, *sources):

        self._point_sources = collections.OrderedDict()
        self._extended_sources = collections.OrderedDict()

        # This is needed to allow the user to access the sources by name as model.source, independently of their
        # type

        self._sources = collections.OrderedDict()

        for source in sources:

            assert source.name not in self._sources, "You cannot use the same name for multiple sources."

            self._sources[source.name] = source

            if isinstance(source, PointSource):

                self._point_sources[source.name] = source

            elif isinstance(source, ExtendedSource):

                self._extended_sources[source.name] = source

        super(Model, self).__init__('source', self._sources)

    @property
    def point_sources(self):
        return self._point_sources

    @property
    def extended_sources(self):
        return self._extended_sources

    @property
    def sources(self):
        return self._sources

    def __repr__(self):

        # Print the name of point sources if there are any

        representation = "Point sources: "

        if self._point_sources:

            representation += ",".join(self._point_sources.keys())

        else:

            representation += "(none)"

        # Print the name of extended sources if there are any

        representation += "\nExtended sources: "

        if self._extended_sources:

            representation += ",".join(self._extended_sources.keys())

        else:

            representation += "(none)"

        return representation

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

        if os.path.exists(output_file) and overwrite==False:

            raise ModelFileExists("The file %s exists already. If you want to overwrite it, use the 'overwrite=True' "
                                  "options as 'model.save(\"%s\", overwrite=True)'. " % (output_file, output_file))

        else:

            # Get the serialization dictionary

            data = self.to_dict()

            # Write it to disk

            try:

                # Get the YAML representation of the data

                representation = my_yaml.dump(data)

                with open(output_file,"w+") as f:

                    # Add a new line at the end of each voice (just for clarity)

                    f.write(representation.replace("\n","\n\n"))

            except IOError:

                raise CannotWriteModel(os.path.dirname(os.path.abspath(output_file)),
                                       "Could not write model file %s. Check your permissions to write or the "
                                       "report on the free space which follows: "% output_file)



