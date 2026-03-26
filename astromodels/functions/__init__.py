# astromodels/functions/__init__.py
from __future__ import annotations

import warnings
from importlib import import_module

# lists in case we get external dependencies also at this level
_exports = {
    # functions_1D
    "Band": (".functions_1D",),
    "Band_Calderone": (".functions_1D",),
    "Band_grbm": (".functions_1D",),
    "Blackbody": (".functions_1D",),
    "Broken_powerlaw": (".functions_1D",),
    "Constant": (".functions_1D",),
    "Cubic": (".functions_1D",),
    "Cutoff_powerlaw": (".functions_1D",),
    "Cutoff_powerlaw_Ep": (".functions_1D",),
    "Cutoff_powerlaw_flux": (".functions_1D",),  # optional (gsl)
    "DiracDelta": (".functions_1D",),
    "DoubleSmoothlyBrokenPowerlaw": (".functions_1D",),
    "Exponential_cutoff": (".functions_1D",),
    "GenericFunction": (".functions_1D",),
    "Inverse_cutoff_powerlaw": (".functions_1D",),
    "Line": (".functions_1D",),
    "Log_parabola": (".functions_1D",),
    "ModifiedBlackbody": (".functions_1D",),
    "NonDissipativePhotosphere": (".functions_1D",),
    "NonDissipativePhotosphere_Deep": (".functions_1D",),
    "PhAbs": (".functions_1D",),
    "Powerlaw": (".functions_1D",),
    "Powerlaw_Eflux": (".functions_1D",),
    "Powerlaw_flux": (".functions_1D",),
    "Quadratic": (".functions_1D",),
    "Quartic": (".functions_1D",),
    "Sin": (".functions_1D",),
    "SmoothlyBrokenPowerLaw": (".functions_1D",),
    "Standard_Rv": (".functions_1D",),
    "StepFunction": (".functions_1D",),
    "StepFunctionUpper": (".functions_1D",),
    "Super_cutoff_powerlaw": (".functions_1D",),
    "TbAbs": (".functions_1D",),
    "WAbs": (".functions_1D",),
    "ZDust": (".functions_1D",),
    "get_polynomial": (".functions_1D",),
    # Optional 1D features (deps handled in functions_1D)
    "APEC": (".functions_1D",),
    "VAPEC": (".functions_1D",),
    "EBLattenuation": (".functions_1D",),
    "Synchrotron": (".functions_1D",),
    # Base function classes
    "Function1D": (".function",),
    "Function2D": (".function",),
    "Function3D": (".function",),
    "FunctionMeta": (".function",),
    # functions_2D
    "Asymm_Gaussian_on_sphere": (".functions_2D",),
    "Disk_on_sphere": (".functions_2D",),
    "Ellipse_on_sphere": (".functions_2D",),
    "Gaussian_on_sphere": (".functions_2D",),
    "Latitude_galactic_diffuse": (".functions_2D",),
    "Power_law_on_sphere": (".functions_2D",),
    "SpatialTemplate_2D": (".functions_2D",),
    # functions_3D
    "Continuous_injection_diffusion": (".functions_3D",),
    "Continuous_injection_diffusion_ellipse": (".functions_3D",),
    "Continuous_injection_diffusion_legacy": (".functions_3D",),
    "GalPropTemplate_3D": (".functions_3D",),
    "Hermes": (".functions_3D",),
    # priors
    "Beta": (".priors",),
    "Cauchy": (".priors",),
    "Cosine_Prior": (".priors",),
    "Exponential": (".priors",),
    "Gamma": (".priors",),
    "Gaussian": (".priors",),
    "Log_normal": (".priors",),
    "Log_uniform_prior": (".priors",),
    "Powerlaw_Prior": (".priors",),
    "Truncated_gaussian": (".priors",),
    "Uniform_prior": (".priors",),
    # template_model
    "MissingDataFile": (".template_model",),
    "TemplateModel": (".template_model",),
    "TemplateModelFactory": (".template_model",),
    "XSPECTableModel": (".template_model",),
    # dark_matter
    "DMFitFunction": (".dark_matter.dm_models",),
    "DMSpectra": (".dark_matter.dm_models",),
    "list_functions": ("astromodels.utils.list_functions",),
    "has_atomdb": (".functions_1D",),
    "has_ebltable": (".functions_1D",),
    "has_naima": (".functions_1D",),
    "has_gsl": (".functions_1D",),
}
_depcrecated = {
    "ModelAssertionViolation": ("astromodels.utils.exceptions",),
}
_exports.update(_depcrecated)

# Public API surface: stable list of names
__all__ = sorted(_exports.keys())


def __getattr__(name: str):
    # Lazy re-export: import the defining module only when the symbol is accessed
    mod_path = _exports.get(name)
    dep_path = _depcrecated.get(name)
    if dep_path:
        warnings.warn(
            f"astromodels.functions.{name} is deprecated; use {dep_path} instead.",
            DeprecationWarning,
            stacklevel=2,
        )
    if mod_path is None:
        # Unknown name
        raise AttributeError(name)
    mod = import_module(mod_path[0], __name__)
    try:
        obj = getattr(mod, name)
    except AttributeError as e:
        raise AttributeError(
            f"{name} is not present in {mod.__name__}. "
            "This may indicate a missing optional dependency or an internal mismatch."
        ) from e
    globals()[name] = obj  # cache binding
    return obj


def __dir__():
    return sorted(__all__)
