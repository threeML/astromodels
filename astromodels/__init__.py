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
                                 SettingOutOfBounds, turn_off_parameter_transforms)
    from .core.polarization import LinearPolarization, StokesPolarization
    from .core.serialization import *
    from .core.spectral_component import SpectralComponent
    from .core.units import get_units
    from .functions import (Asymm_Gaussian_on_sphere, Band, Band_Calderone,
                            Band_grbm, Blackbody, Broken_powerlaw, Cauchy,
                            Constant, Continuous_injection_diffusion,
                            Continuous_injection_diffusion_ellipse,
                            Continuous_injection_diffusion_legacy,
                            Cosine_Prior, Cubic, Cutoff_powerlaw,
                            Cutoff_powerlaw_Ep, DiracDelta, Disk_on_sphere,
                            DMFitFunction, DMSpectra, Ellipse_on_sphere,
                            Exponential_cutoff, Function1D, Function2D,
                            Function3D, FunctionMeta, GalPropTemplate_3D,
                            Gaussian, Gaussian_on_sphere,
                            Inverse_cutoff_powerlaw, Latitude_galactic_diffuse,
                            Line, Log_normal, Log_parabola, Log_uniform_prior,
                            MissingDataFile, ModelAssertionViolation,
                            ModifiedBlackbody, NonDissipativePhotosphere,
                            NonDissipativePhotosphere_Deep, PhAbs,
                            Power_law_on_sphere, Powerlaw, Powerlaw_Eflux,
                            Powerlaw_flux, Quadratic, Quartic, Sin,
                            SmoothlyBrokenPowerLaw, SpatialTemplate_2D,
                            Standard_Rv, StepFunction, StepFunctionUpper,
                            Super_cutoff_powerlaw, TbAbs, TemplateModel,
                            TemplateModelFactory, Truncated_gaussian,
                            Uniform_prior, WAbs, XSPECTableModel, ZDust,
                            get_polynomial, has_ebltable, has_gsl, has_naima, has_atomdb)

    if has_ebltable:

        from .functions import EBLattenuation

    if has_gsl:

        from .functions import Cutoff_powerlaw_flux

    if has_naima:

        from .functions import Synchrotron

    if has_atomdb:

        from .functions import APEC, VAPEC
        
    from .functions.function import get_function_class, list_functions
    from .sources import ExtendedSource, ParticleSource, PointSource


    astromodels_units = get_units()
    from astromodels.utils.logging import setup_logger, update_logging_level, silence_warnings, activate_warnings

import astropy.units as u

log = setup_logger(__name__)

__version__ = get_versions()['version']
del get_versions
