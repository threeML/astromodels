# Import the version
from version import __version__

from .sources.point_source import PointSource
from .sources.extended_source import ExtendedSource
from .sources.particle_source import ParticleSource
from parameter import Parameter, IndependentVariable, SettingOutOfBounds
from .functions.functions import *
from .functions.functions_2D import *
from .functions.functions_3D import *
from .functions.function import list_functions, get_function_class
from model import Model
from spectral_component import SpectralComponent
from model_parser import load_model
from units import get_units
from .xspec.factory import setup_xspec_models, has_xspec
import astropy.units as u

astromodels_units = get_units()

if has_xspec:

    from .xspec.factory import *
    from .xspec.xspec_settings import *