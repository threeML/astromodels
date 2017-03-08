# Import the version
# from version import __version__
#
#

import os

if os.environ.get("ASTROMODELS_DEBUG", None) is None:

    from .sources.point_source import PointSource
    from .sources.extended_source import ExtendedSource
    from .sources.particle_source import ParticleSource
    from .core.parameter import Parameter, IndependentVariable, SettingOutOfBounds
    from .functions.functions import *
    from .functions.functions_2D import *
    from .functions.functions_3D import *
    from .functions.template_model import *
    from .functions. dark_matter.dm_models import *
    from .functions.function import list_functions, get_function_class
    from .core.model import Model
    from .core.parameter import Parameter
    from .core.spectral_component import SpectralComponent
    from .core.model_parser import load_model, clone_model
    from .core.units import get_units
    from .core.memoization import use_astromodels_memoization

    astromodels_units = get_units()


import astropy.units as u

# if has_xspec:
#
#     from .xspec.factory import *
#     from .xspec.xspec_settings import *