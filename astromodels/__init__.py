# astromodels/__init__.py
from __future__ import annotations

import importlib
import warnings

from ._version import get_versions

__version__ = get_versions()["version"]
del get_versions

# Map top-level public names to (module, attribute)
# These include the legacy function/prior names so that __all__ matches historical
# behavior.
_functions = {
    # 1D functions
    "Blackbody": ("astromodels.functions.functions_1D.blackbody", "Blackbody"),
    "ModifiedBlackbody": (
        "astromodels.functions.functions_1D.blackbody",
        "ModifiedBlackbody",
    ),
    "NonDissipativePhotosphere": (
        "astromodels.functions.functions_1D.blackbody",
        "NonDissipativePhotosphere",
    ),
    "NonDissipativePhotosphere_Deep": (
        "astromodels.functions.functions_1D.blackbody",
        "NonDissipativePhotosphere_Deep",
    ),
    "GenericFunction": (
        "astromodels.functions.functions_1D.functions",
        "GenericFunction",
    ),
    "StepFunction": ("astromodels.functions.functions_1D.functions", "StepFunction"),
    "StepFunctionUpper": (
        "astromodels.functions.functions_1D.functions",
        "StepFunctionUpper",
    ),
    "Sin": ("astromodels.functions.functions_1D.functions", "Sin"),
    "DiracDelta": ("astromodels.functions.functions_1D.functions", "DiracDelta"),
    "Log_parabola": ("astromodels.functions.functions_1D.functions", "Log_parabola"),
    "Exponential_cutoff": (
        "astromodels.functions.functions_1D.functions",
        "Exponential_cutoff",
    ),
    "PhAbs": ("astromodels.functions.functions_1D.absorption", "PhAbs"),
    "TbAbs": ("astromodels.functions.functions_1D.absorption", "TbAbs"),
    "WAbs": ("astromodels.functions.functions_1D.absorption", "WAbs"),
    "ZDust": ("astromodels.functions.functions_1D.extinction", "ZDust"),
    "Constant": ("astromodels.functions.functions_1D.polynomials", "Constant"),
    "Line": ("astromodels.functions.functions_1D.polynomials", "Line"),
    "Quadratic": ("astromodels.functions.functions_1D.polynomials", "Quadratic"),
    "Cubic": ("astromodels.functions.functions_1D.polynomials", "Cubic"),
    "Quartic": ("astromodels.functions.functions_1D.polynomials", "Quartic"),
    "Powerlaw": ("astromodels.functions.functions_1D.powerlaws", "Powerlaw"),
    "Powerlaw_flux": ("astromodels.functions.functions_1D.powerlaws", "Powerlaw_flux"),
    "Powerlaw_Eflux": (
        "astromodels.functions.functions_1D.powerlaws",
        "Powerlaw_Eflux",
    ),
    "Cutoff_powerlaw": (
        "astromodels.functions.functions_1D.powerlaws",
        "Cutoff_powerlaw",
    ),
    "Cutoff_powerlaw_Ep": (
        "astromodels.functions.functions_1D.powerlaws",
        "Cutoff_powerlaw_Ep",
    ),
    "Inverse_cutoff_powerlaw": (
        "astromodels.functions.functions_1D.powerlaws",
        "Inverse_cutoff_powerlaw",
    ),
    "Super_cutoff_powerlaw": (
        "astromodels.functions.functions_1D.powerlaws",
        "Super_cutoff_powerlaw",
    ),
    "SmoothlyBrokenPowerLaw": (
        "astromodels.functions.functions_1D.powerlaws",
        "SmoothlyBrokenPowerLaw",
    ),
    "Broken_powerlaw": (
        "astromodels.functions.functions_1D.powerlaws",
        "Broken_powerlaw",
    ),
    "Band": ("astromodels.functions.functions_1D.powerlaws", "Band"),
    "Band_grbm": ("astromodels.functions.functions_1D.powerlaws", "Band_grbm"),
    "Band_Calderone": (
        "astromodels.functions.functions_1D.powerlaws",
        "Band_Calderone",
    ),
    "DoubleSmoothlyBrokenPowerlaw": (
        "astromodels.functions.functions_1D.powerlaws",
        "DoubleSmoothlyBrokenPowerlaw",
    ),
    # 2D functions
    "Latitude_galactic_diffuse": (
        "astromodels.functions.functions_2D",
        "Latitude_galactic_diffuse",
    ),
    "Gaussian_on_sphere": ("astromodels.functions.functions_2D", "Gaussian_on_sphere"),
    "Asymm_Gaussian_on_sphere": (
        "astromodels.functions.functions_2D",
        "Asymm_Gaussian_on_sphere",
    ),
    "Disk_on_sphere": ("astromodels.functions.functions_2D", "Disk_on_sphere"),
    "Ellipse_on_sphere": ("astromodels.functions.functions_2D", "Ellipse_on_sphere"),
    "SpatialTemplate_2D": ("astromodels.functions.functions_2D", "SpatialTemplate_2D"),
    "Power_law_on_sphere": (
        "astromodels.functions.functions_2D",
        "Power_law_on_sphere",
    ),
    # 3D functions
    "Continuous_injection_diffusion_ellipse": (
        "astromodels.functions.functions_3D",
        "Continuous_injection_diffusion_ellipse",
    ),
    "Continuous_injection_diffusion": (
        "astromodels.functions.functions_3D",
        "Continuous_injection_diffusion",
    ),
    "Continuous_injection_diffusion_legacy": (
        "astromodels.functions.functions_3D",
        "Continuous_injection_diffusion_legacy",
    ),
    "Hermes": ("astromodels.functions.functions_3D", "Hermes"),
    # Priors
    "Gaussian": ("astromodels.functions.priors", "Gaussian"),
    "Truncated_gaussian": ("astromodels.functions.priors", "Truncated_gaussian"),
    "Cauchy": ("astromodels.functions.priors", "Cauchy"),
    "Cosine_Prior": ("astromodels.functions.priors", "Cosine_Prior"),
    "Log_normal": ("astromodels.functions.priors", "Log_normal"),
    "Uniform_prior": ("astromodels.functions.priors", "Uniform_prior"),
    "Log_uniform_prior": ("astromodels.functions.priors", "Log_uniform_prior"),
    "Beta": ("astromodels.functions.priors", "Beta"),
    "Gamma": ("astromodels.functions.priors", "Gamma"),
    "Exponential": ("astromodels.functions.priors", "Exponential"),
    "Powerlaw_Prior": ("astromodels.functions.priors", "Powerlaw_Prior"),
    # Templates
    "TemplateModel": ("astromodels.functions.template_model", "TemplateModel"),
    "TemplateModelFactory": (
        "astromodels.functions.template_model",
        "TemplateModelFactory",
    ),
    "XSPECTableModel": ("astromodels.functions.template_model", "XSPECTableModel"),
    "MissingDataFile": ("astromodels.functions.template_model", "MissingDataFile"),
    # Dark matter
    "DMFitFunction": ("astromodels.functions.dark_matter.dm_models", "DMFitFunction"),
    "DMSpectra": ("astromodels.functions.dark_matter.dm_models", "DMSpectra"),
    "show_configuration": ("astromodels.utils.configuration", "show_configuration"),
    "functions": ("astromodels.functions", "functions"),
    "list_functions": ("astromodels.functions.function", "list_functions"),
}

_public = {
    # core
    "Model": ("astromodels.core.model", "Model"),
    "clone_model": ("astromodels.core.model_parser", "clone_model"),
    "load_model": ("astromodels.core.model_parser", "load_model"),
    "IndependentVariable": ("astromodels.core.parameter", "IndependentVariable"),
    "Parameter": ("astromodels.core.parameter", "Parameter"),
    "turn_off_parameter_transforms": (
        "astromodels.core.parameter",
        "turn_off_parameter_transforms",
    ),
    "LinearPolarization": ("astromodels.core.polarization", "LinearPolarization"),
    "StokesPolarization": ("astromodels.core.polarization", "StokesPolarization"),
    "serialize_model": ("astromodels.core.serialization", "serialize_model"),
    "unserialize_model": ("astromodels.core.serialization", "unserialize_model"),
    "SpectralComponent": ("astromodels.core.spectral_component", "SpectralComponent"),
    "use_astromodels_memoization": (
        "astromodels.core.memoization",
        "use_astromodels_memoization",
    ),
    "ModelAssertionViolation": (
        "astromodels.utils.exceptions",
        "ModelAssertionViolation",
    ),
    "SettingOutOfBounds": ("astromodels.core.parameter", "SettingOutOfBounds"),
    # sources
    "ExtendedSource": ("astromodels.sources.extended_source", "ExtendedSource"),
    "ParticleSource": ("astromodels.sources.particle_source", "ParticleSource"),
    "PointSource": ("astromodels.sources.point_source", "PointSource"),
}

# Merge everything we want to expose at top level
_public.update(_functions)

# Which top-level names are deprecated (emit warnings when accessed)
DEPRECATED_TOPLEVEL = set(_functions.keys())

# Export everything (historic behavior), plus convenience names
__all__ = sorted(set(_public.keys()) | {"__version__", "astromodels_config"})


def __getattr__(name: str):
    # Lazy re-exports
    if name == "astromodels_config":
        mod = importlib.import_module("astromodels.utils.configuration")
        val = getattr(mod, "astromodels_config")
        globals()[name] = val
        return val

    try:
        mod_name, attr = _public[name]
    except KeyError:
        raise AttributeError(
            f"module 'astromodels' has no attribute {name!r}"
        ) from None

    # Emit deprecation for legacy top-level function/prior/template names
    if name in DEPRECATED_TOPLEVEL:
        warnings.warn(
            f"Top-level access 'astromodels.{name}' is deprecated; "
            f"use 'from {mod_name} import {name}' instead.",
            category=DeprecationWarning,
            stacklevel=2,
        )

    try:
        mod = importlib.import_module(mod_name)
        val = getattr(mod, attr)
    except ImportError as e:
        # Surface clearer message for optional dependencies or missing submodules
        raise ImportError(
            f"Cannot access 'astromodels.{name}' because '{mod_name}' "
            f"could not be imported. This feature may require optional dependencies. "
            f"Original error: {e}"
        ) from e
    except AttributeError as e:
        # Defensive: module imported but symbol missing (internal mismatch)
        raise AttributeError(
            f"'astromodels.{name}' could not be resolved from '{mod_name}'."
        ) from e

    globals()[name] = val  # cache
    return val


def __dir__():
    # List everything we export (historic behavior)
    return sorted(__all__)
