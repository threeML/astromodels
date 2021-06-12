import collections

import astropy.units as u
import numpy as np

from astromodels.core.spectral_component import SpectralComponent
from astromodels.core.tree import Node
from astromodels.core.units import get_units
from astromodels.functions import Constant
from astromodels.sources.source import Source, SourceType
from astromodels.utils.pretty_list import dict_to_list
from astromodels.utils.logging import setup_logger

log = setup_logger(__name__)

class ExtendedSource(Source, Node):
    def __init__(self,
                 source_name,
                 spatial_shape,
                 spectral_shape=None,
                 components=None):

        # Check that we have all the required information
        # and set the units

        current_u = get_units()

        if spatial_shape.n_dim == 2:

            # Now gather the component(s)

            # We need either a single component, or a list of components, but not both
            # (that's the ^ symbol)

            assert (spectral_shape is not None) ^ (components is not None), "You have to provide either a single " \
                                                                            "component, or a list of components " \
                                                                            "(but not both)."

            # If the user specified only one component, make a list of one element with a default name ("main")

            if spectral_shape is not None:

                components = [SpectralComponent("main", spectral_shape)]

            # Components in this case have energy as x and differential flux as y

            diff_flux_units = (current_u.energy * current_u.area *
                               current_u.time)**(-1)

            # Now set the units of the components
            for component in components:

                component.shape.set_units(current_u.energy, diff_flux_units)

            # Set the units of the brightness
            spatial_shape.set_units(current_u.angle, current_u.angle,
                                    current_u.angle**(-2))

        elif spatial_shape.n_dim == 3:

            # If there is no spectral component then assume that the input is a template, which will provide the
            # spectrum by itself. We just use a renormalization (a bias)

            if spectral_shape is None and components is None:

                # This is a template. Add a component which is just a renormalization

                spectral_shape = Constant()
                components = [SpectralComponent("main", spectral_shape)]

                # set the units
                diff_flux_units = (current_u.energy * current_u.area *
                                   current_u.time * current_u.angle**2)**(-1)
                spatial_shape.set_units(current_u.angle, current_u.angle,
                                        current_u.energy, diff_flux_units)

            else:

                # the spectral shape has been given, so this is a case where the spatial template gives an
                # energy-dependent shape and the spectral components give the spectrum

                if not ((spectral_shape is not None) ^ (components is not None)):

                    log.error("You can provide either a single "
                              "component, or a list of components "
                              "(but not both).")

                    raise AssertionError()
                    
                if spectral_shape is not None:

                    components = [SpectralComponent("main", spectral_shape)]

                # Assign units
                diff_flux_units = (current_u.energy * current_u.area * current_u.time) ** (-1)

                # Now set the units of the components
                for component in components:
                    component.shape.set_units(current_u.energy, diff_flux_units)

                # Set the unit of the spatial template
                spatial_shape.set_units(current_u.angle, current_u.angle, current_u.energy, current_u.angle**(-2))

        else:

            log.error("The spatial shape must have either 2 or 3 dimensions.")

            raise RuntimeError()

        # Here we have a list of components

        Source.__init__(self, components, SourceType.EXTENDED_SOURCE)

        # A source is also a Node in the tree

        Node.__init__(self, source_name)

        # Add the spatial shape as a child node, with an explicit name
        self._spatial_shape = spatial_shape
        self._add_child(self._spatial_shape)

        # Add the same node also with the name of the function
        #self._add_child(self._shape, self._shape.__name__)

        # Add a node called 'spectrum'

        spectrum_node = Node('spectrum')
        spectrum_node._add_children(list(self._components.values()))

        self._add_child(spectrum_node)

    @property
    def spatial_shape(self):
        """
        A generic name for the spatial shape.

        :return: the spatial shape instance
        """

        return self._spatial_shape

    def get_spatially_integrated_flux( self, energies):
    
        """
        Returns total flux of source at the given energy
        :param energies: energies (array or float)
        :return: differential flux at given energy
        """
       
        if not isinstance(energies, np.ndarray):
            energies = np.array(energies, ndmin=1)

        # Get the differential flux from the spectral components

        results = [self.spatial_shape.get_total_spatial_integral(energies) * component.shape(energies) for component in self.components.values()]

        if isinstance(energies, u.Quantity):

            # Slow version with units

            # We need to sum like this (slower) because using np.sum will not preserve the units
            # (thanks astropy.units)

            differential_flux = sum(results)

        else:

            # Fast version without units, where x is supposed to be in the same units as currently defined in
            # units.get_units()

            differential_flux = np.sum(results, 0)

        return differential_flux


    def __call__(self, lon, lat, energies):
        """
        Returns brightness of source at the given position and energy
        :param lon: longitude (array or float)
        :param lat: latitude (array or float)
        :param energies: energies (array or float)
        :return: differential flux at given position and energy
        """

        assert type(lat) == type(lon) and type(lon) == type(energies), "Type mismatch in input of call"

        if not isinstance(lat, np.ndarray):

            lat = np.array(lat, ndmin=1)
            lon = np.array(lon, ndmin=1)
            energies = np.array(energies, ndmin=1)

        # Get the differential flux from the spectral components

        results = [component.shape(energies) for component in list(self.components.values())]

        if isinstance(energies, u.Quantity):

            # Slow version with units

            # We need to sum like this (slower) because using np.sum will not preserve the units
            # (thanks astropy.units)

            differential_flux = sum(results)

        else:

            # Fast version without units, where x is supposed to be in the same units as currently defined in
            # units.get_units()

            differential_flux = np.sum(results, 0)

        # Get brightness from spatial model

        if self._spatial_shape.n_dim == 2:

            brightness = self._spatial_shape(lon, lat)

            # In this case the spectrum is the same everywhere
            n_points = lat.shape[0]
            n_energies = differential_flux.shape[0]

            # The following is a little obscure, but it is 6x faster than doing a for loop

            cube = np.repeat(differential_flux, n_points).reshape(n_energies, n_points).T
            result = (cube.T * brightness).T

        else:

            result = self._spatial_shape(lon, lat, energies) * differential_flux

        # Do not clip the output, otherwise it will not be possible to use ext. sources
        # with negative fluxes

        return np.squeeze(result)

    @property
    def has_free_parameters(self):
        """
        Returns True or False whether there is any parameter in this source

        :return:
        """

        for component in list(self._components.values()):

            for par in list(component.shape.parameters.values()):

                if par.free:

                    return True

        for par in list(self.spatial_shape.parameters.values()):

            if par.free:

                return True

        return False

    @property
    def free_parameters(self):
        """
        Returns a dictionary of free parameters for this source
        We use the parameter path as the key because it's 
        guaranteed to be unique, unlike the parameter name.

        :return:
        """
        free_parameters = collections.OrderedDict()

        for component in list(self._components.values()):

            for par in list(component.shape.parameters.values()):

                if par.free:

                    free_parameters[par.path] = par

        for par in list(self.spatial_shape.parameters.values()):

            if par.free:

                free_parameters[par.path] = par

        return free_parameters

    @property
    def parameters(self):
        """
        Returns a dictionary of all parameters for this source.
        We use the parameter path as the key because it's 
        guaranteed to be unique, unlike the parameter name.

        :return:
        """
        all_parameters = collections.OrderedDict()

        for component in list(self._components.values()):

            for par in list(component.shape.parameters.values()):

                all_parameters[par.path] = par

        for par in list(self.spatial_shape.parameters.values()):

            all_parameters[par.path] = par

        return all_parameters

    def _repr__base(self, rich_output=False):
        """
        Representation of the object

        :param rich_output: if True, generates HTML, otherwise text
        :return: the representation
        """

        # Make a dictionary which will then be transformed in a list

        repr_dict = collections.OrderedDict()

        key = '%s (extended source)' % self.name

        repr_dict[key] = collections.OrderedDict()
        repr_dict[key]['shape'] = self._spatial_shape.to_dict(minimal=True)
        repr_dict[key]['spectrum'] = collections.OrderedDict()

        for component_name, component in list(self.components.items()):
            repr_dict[key]['spectrum'][component_name] = component.to_dict(minimal=True)

        return dict_to_list(repr_dict, rich_output)

    def get_boundaries(self):
        """
        Returns the boundaries for this extended source

        :return: a tuple of tuples ((min. lon, max. lon), (min lat, max lat))
        """
        return self._spatial_shape.get_boundaries()
