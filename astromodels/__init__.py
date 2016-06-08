from .sources.point_source import PointSource
from .sources.extended_source import ExtendedSource
from .sources.particle_source import ParticleSource
from parameter import Parameter, IndependentVariable
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

# This will either work or issue a warning if XSpec is not available

new_functions = setup_xspec_models()

# Now import the new classes in the local namespace (if any)
for function_name in new_functions:

    locals()[function_name] = get_function_class(function_name)

if has_xspec:

    from .xspec.xspec_settings import *