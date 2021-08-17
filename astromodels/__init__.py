from __future__ import absolute_import

import os

from ._version import get_versions

# Import the version

#
#

if os.environ.get("ASTROMODELS_DEBUG", None) is None:

    from .utils.configuration import astromodels_config
    from .core.memoization import use_astromodels_memoization
    from .core.model import Model
    from .core.model_parser import clone_model, load_model
    from .core.parameter import (IndependentVariable, Parameter,
                                 SettingOutOfBounds)
    from .core.polarization import LinearPolarization, StokesPolarization
    from .core.serialization import *
    from .core.spectral_component import SpectralComponent
    from .core.units import get_units
    from .functions import (
        Band, Band_Calderone, Band_grbm, Broken_powerlaw, Cutoff_powerlaw,
        Inverse_cutoff_powerlaw, Powerlaw, Powerlaw_Eflux, Powerlaw_flux,
        SmoothlyBrokenPowerLaw, Super_cutoff_powerlaw, Constant, Cubic,
        DiracDelta, Exponential_cutoff, Line, Quadratic, Sin, StepFunction,
        StepFunctionUpper, PhAbs, TbAbs, WAbs, Asymm_Gaussian_on_sphere,
        Disk_on_sphere, Ellipse_on_sphere, Gaussian_on_sphere,
        Latitude_galactic_diffuse, Power_law_on_sphere, SpatialTemplate_2D,
        Continuous_injection_diffusion, Continuous_injection_diffusion_ellipse,
        Continuous_injection_diffusion_legacy, GalPropTemplate_3D, DMSpectra,
        DMFitFunction, Cauchy, Cosine_Prior, Gaussian, Log_normal, Quartic,
        get_polynomial, Log_uniform_prior, Truncated_gaussian, Uniform_prior,
        TemplateModel, TemplateModelFactory, XSPECTableModel, MissingDataFile,
        Log_parabola, Blackbody, Function1D, Function2D, Function3D,
        FunctionMeta, ModelAssertionViolation, has_naima, has_gsl,
        has_ebltable, has_atomdb)

    if has_ebltable:

        from .functions import EBLattenuation

    if has_gsl:

        from .functions import Cutoff_powerlaw_flux

    if has_naima:

        from .functions import Synchrotron

    if has_atomdb:

        from .functions import APEC, VAPEC
        
    from .functions.function import get_function_class, list_functions
    from .sources import ExtendedSource, PointSource, ParticleSource


    astromodels_units = get_units()
    from astromodels.utils.logging import setup_logger, update_logging_level, silence_warnings, activate_warnings

import astropy.units as u

log = setup_logger(__name__)

__version__ = get_versions()['version']
del get_versions
