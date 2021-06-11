from .functions_1D import (
    Blackbody, Constant, Cubic, DiracDelta, Exponential_cutoff, Line,
    Log_parabola, Quadratic, Sin, StepFunction, StepFunctionUpper,
    has_ebltable, has_gsl, has_naima, Band, Band_Calderone, Band_grbm,
    Broken_powerlaw, Cutoff_powerlaw, Cutoff_powerlaw_Ep, Inverse_cutoff_powerlaw, Powerlaw,
    Powerlaw_Eflux, Powerlaw_flux, SmoothlyBrokenPowerLaw, Quartic,
    get_polynomial, Super_cutoff_powerlaw, PhAbs, TbAbs, WAbs, has_atomdb)

if has_naima:
    from .functions_1D import Synchrotron

if has_gsl:

    from .functions_1D import Cutoff_powerlaw_flux

if has_ebltable:
    from .functions_1D import EBLattenuation

if has_atomdb:

    from .apec import APEC, VAPEC

from .dark_matter.dm_models import DMFitFunction, DMSpectra

from .function import (Function1D, Function2D, Function3D, FunctionMeta,
                       ModelAssertionViolation)

from .functions_2D import (Asymm_Gaussian_on_sphere, Disk_on_sphere,
                           Ellipse_on_sphere, Gaussian_on_sphere,
                           Latitude_galactic_diffuse, Power_law_on_sphere,
                           SpatialTemplate_2D)
from .functions_3D import (Continuous_injection_diffusion,
                           Continuous_injection_diffusion_ellipse,
                           Continuous_injection_diffusion_legacy,
                           GalPropTemplate_3D)
from .priors import (Cauchy, Cosine_Prior, Gaussian, Log_normal,
                     Log_uniform_prior, Truncated_gaussian, Uniform_prior)
from .template_model import (MissingDataFile, TemplateModel,
                             TemplateModelFactory, XSPECTableModel)

__all__ = [
    "Band", "Band_Calderone", "Band_grbm", "Broken_powerlaw",
    "Cutoff_powerlaw", "Cutoff_powerlaw_Ep", "Inverse_cutoff_powerlaw", "Powerlaw", "Powerlaw_Eflux",
    "Powerlaw_flux", "SmoothlyBrokenPowerLaw", "Super_cutoff_powerlaw",
    "Constant", "Cubic", "DiracDelta", "Exponential_cutoff", "Line",
    "Quadratic", "Sin", "StepFunction", "StepFunctionUpper", "PhAbs", "TbAbs",
    "WAbs", "Asymm_Gaussian_on_sphere", "Disk_on_sphere", "Ellipse_on_sphere",
    "Gaussian_on_sphere", "Latitude_galactic_diffuse", "Power_law_on_sphere",
    "SpatialTemplate_2D", "Continuous_injection_diffusion",
    "Continuous_injection_diffusion_ellipse",
    "Continuous_injection_diffusion_legacy", "GalPropTemplate_3D", "DMSpectra",
    "DMFitFunction", "Cauchy", "Cosine_Prior", "Gaussian", "Log_normal",
    "Log_uniform_prior", "Truncated_gaussian", "Uniform_prior",
    "TemplateModel", "TemplateModelFactory", "XSPECTableModel",
    "MissingDataFile", "Log_parabola", "Blackbody", "Function1D", "Function2D",
    "Function3D", "FunctionMeta", "ModelAssertionViolation", "Quartic", "get_polynomial"
]

if has_atomdb:
    __all__.extend(["APEC", "VAPEC"])

if has_gsl:

    __all__.extend(["Cutoff_powerlaw_flux"])

if has_naima:

    __all__.extend(["Synchrotron"])

if has_ebltable:

    __all__.extend(["EBLattenuation"])
