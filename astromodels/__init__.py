from __future__ import absolute_import

import os

from ._version import get_versions

# Import the version

#
#


if os.environ.get("ASTROMODELS_DEBUG", None) is None:

    from .core.memoization import use_astromodels_memoization
    from .core.model import Model
    from .core.model_parser import clone_model, load_model
    from .core.parameter import (IndependentVariable, Parameter,
                                 SettingOutOfBounds)
    from .core.polarization import LinearPolarization, StokesPolarization
    from .core.serialization import *
    from .core.spectral_component import SpectralComponent
    from .core.units import get_units
    from .functions import *
    from .functions.function import get_function_class, list_functions
    from .sources.extended_source import ExtendedSource
    from .sources.particle_source import ParticleSource
    from .sources.point_source import PointSource

    astromodels_units = get_units()
    from astromodels.utils.logging import setup_logger, update_logging_level, silence_warnings, activate_warnings

import astropy.units as u

log = setup_logger(__name__)

__version__ = get_versions()['version']
del get_versions
