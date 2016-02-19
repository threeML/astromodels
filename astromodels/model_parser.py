__author__ = 'giacomov'

__doc__ = """
=============
Description
=============

Models in astromodels are defined with a language as closest as possible to the natural language. In the definition of
the model, which happens only once in any given analysis, the syntax is very explicit and analytical, while the use of
the Model class is geared towards usability.

This is an example of a simple likelihood model, with one point source having two spectral components, both with
a powerlaw shape but the first with a linear polarization, the second with a circular polarization::

  test_source :

      point source:

          position:

              RA : {value: 23.5}
              Dec : {value: -45.4}

          spectrum :

              main :

                  shape:

                      powerlaw :

                          logK : {value: 1}
                          index : {value: -2}
                          piv : {value: 100, unit: 'keV'}

                  polarization :

                      linear :

                         degree : {value: 0.2}
                         angle : {value: 45}

              other :

                  shape:

                      powerlaw :

                          logK : {value: 1}
                          index : {value: -2}
                          piv : {value: 200}

                  polarization :

                      stokes:

                          I : {value : 1}
                          Q : {value : 0}
                          U : {value : 0}
                          V : {value : 1}

This model can be read like this:

  > mp = model_parser.ModelParser('test.yml')
  > mod = m.get_model()

Now mod is a instance of a Model, and parameters can be accessed like::

  > ra = mod.test_source.position.RA
  > dec = mod.test_source.position.Dec
  > logK = mod.test_source.main.shape.logK
  > I = mod.test_source.other.polarization.I

Although perhaps a more natural way would have been "mod.test_source.spectrum.main.shape.powerlaw.logK, this would have
resulted in a very long sequence. Thus, all redundant expressions have been removed.

"""

from astromodels import sky_direction
from astromodels.functions import function
from astromodels import spectral_component
from astromodels import polarization
from astromodels.sources import point_source
from astromodels import parameter
from astromodels import model
from astromodels.my_yaml import my_yaml


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


class ModelParser(object):

    def __init__(self, model_file):

        # Read model file and deserialize into a dictionary

        try:

            with open(model_file) as f:

                self._model_dict = my_yaml.load(f)

        except IOError:

            raise ModelIOError("File %s cannot be read. Check path and permissions for current user." % model_file)

        except my_yaml.YAMLError:

            raise ModelYAMLError("Could not parse file %s. Check your syntax." % model_file)

        self._parse()

    def _parse(self):

        # Traverse the dictionary and create all the needed classes

        # The first level is the source level

        self._sources = []

        for source_name, source_definition in self._model_dict.iteritems():

            this_parser = SourceParser(source_name, source_definition)

            self._sources.append(this_parser.get_source())

    def get_model(self):

        return model.Model(*self._sources)


class SourceParser(object):

    def __init__(self, source_name, source_definition):

        assert len(source_definition.values())==1, "Source %s cannot have multiple types" % source_name

        self._source_name = source_name

        # Point source or extended source?

        source_type = source_definition.keys()[0]

        if source_type == 'point source':

            self._parsed_source = self._parse_point_source(source_definition['point source'])

        elif source_type == 'extended source':

            self._parsed_source = self._parse_extended_source(source_definition['extended source'])

        else:

            raise ModelSyntaxError("Don't recognize source type '%s' for source %s. "
                                   "Valid types are 'point source' or 'extended source'." %(source_type, source_name) )

    def get_source(self):

        return self._parsed_source

    def _parse_point_source(self, pts_source_definition):

        # Parse the positional information

        try:

            position_definition = pts_source_definition['position']

        except KeyError:

            raise ModelSyntaxError("Point source %s is missing the 'position' attribute" % self._source_name)

        this_sky_direction = self._parse_sky_direction(position_definition)

        # Parse the spectral information

        try:

            spectrum = pts_source_definition['spectrum']

        except KeyError:

            raise ModelSyntaxError("Point source %s is missing the 'spectrum' attribute" % self._source_name)

        components = []

        for component_name, component_definition in pts_source_definition['spectrum'].iteritems():

            this_component = self._parse_spectral_component(component_name, component_definition)

            components.append(this_component)

        this_point_source = point_source.PointSource(self._source_name, sky_position=this_sky_direction,
                                                     components=components)

        return this_point_source

    def _parse_sky_direction(self, sky_direction_definition):

        # Instance the SkyDirection class using the coordinates provided

        coordinates = {}

        if 'ra' in sky_direction_definition and 'dec' in sky_direction_definition:

            ra = parameter.Parameter('RA',sky_direction_definition['ra']['value'])
            ra.set_bounds(0, 360)
            ra.fix = True

            dec = parameter.Parameter('Dec',sky_direction_definition['dec']['value'])
            dec.set_bounds(-90, 90)
            dec.fix = True

            coordinates['ra'] = ra
            coordinates['dec'] = dec

        elif 'l' in sky_direction_definition and 'b' in sky_direction_definition:

            l = parameter.Parameter('l',sky_direction_definition['l']['value'])
            l.set_bounds(0, 360)
            l.fix = True

            b = parameter.Parameter('b',sky_direction_definition['b']['value'])
            b.set_bounds(-90, 90)
            b.fix = True

            coordinates['l'] = l
            coordinates['b'] = b

        else:

            raise ModelSyntaxError("Position specification for source %s has an invalid coordinate pair. "
                                   " You need to specify either 'ra' and 'dec', or 'l' and 'b'."
                                   % self._source_name)

        # Check if there is a equinox specification

        if 'equinox' in sky_direction_definition:

            coordinates['equinox'] = sky_direction_definition['equinox']

        try:

            this_sky_direction = sky_direction.SkyDirection(**coordinates)

        except sky_direction.WrongCoordinatePair:

            raise ModelSyntaxError("Position specification for source %s has an invalid coordinate pair"
                                   % self._source_name)

        return this_sky_direction

    def _parse_spectral_component(self, component_name, component_definition):

        # Parse the shape definition

        try:

            shape_definition = component_definition['shape']

        except KeyError:

            raise ModelSyntaxError("The component %s of source %s is missing the 'shape' attribute"
                                   % (component_name, self._source_name))

        shape = self._parse_shape_definition(component_name, shape_definition)

        # this_polarization = polarization.Polarization() #TODO

        this_polarization = None

        this_spectral_component = spectral_component.SpectralComponent(component_name, shape, this_polarization)

        return this_spectral_component

    def _parse_shape_definition(self, component_name, shape_definition):

        assert len(shape_definition.values())==1, "Source %s has more than one shape" % self._source_name

        # Get the name of the function

        function_name = shape_definition.keys()[0]

        try:

            this_function = function.get_function(function_name)

        except AttributeError:

            raise ModelSyntaxError("Function %s, specified as shape for component %s of source %s, is not a "
                                   "known 1d function" %(function_name, component_name, self._source_name))

        # Instance the function and set its parameters

        function_instance = this_function()

        parameters_definition = shape_definition[function_name]

        # Loop over the parameters of the function instance, instead of the specification,
        # so we can understand if there are parameters missing from the specification

        for parameter_name in function_instance.parameters.keys():

            try:

                this_definition = parameters_definition[parameter_name]

            except KeyError:

                raise ModelSyntaxError("Function %s, specified as shape for component %s of source %s, lacks "
                                       "the definition for parameter %s"
                                       %(function_name, component_name, self._source_name, parameter_name))

            # Update the parameter. Note that the order is important, because trying to set the value before the
            # minimum and maximum could result in a error.

            # All these specifications are optional. If they are not present, then the default value
            # already contained in the instance of the function will be used

            if 'min' in this_definition:

                function_instance.parameters[parameter_name].min_value = this_definition['min']

            if 'max' in this_definition:

                function_instance.parameters[parameter_name].max_value = this_definition['max']

            if 'delta' in this_definition:

                function_instance.parameters[parameter_name].delta = this_definition['delta']

            if 'fix' in this_definition:

                function_instance.parameters[parameter_name].fix = this_definition['fix']

            if 'free' in this_definition:

                function_instance.parameters[parameter_name].free = this_definition['free']

            if 'unit' in this_definition:

                function_instance.parameters[parameter_name].unit = this_definition['unit']

            # Now set the value, which must be present

            try:

                function_instance.parameters[parameter_name].value = this_definition['value']

            except KeyError:

                raise ModelSyntaxError("The parameter %s in function %s, specified as shape for component %s "
                                       "of source %s, lacks a 'value' attribute"
                                       % (parameter_name, function_name, component_name, self._source_name))
        return function_instance

    def _parse_extended_source(self, ext_source_definition):

        return 0