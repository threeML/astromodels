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
_func_names = [
    "Blackbody",
    "ModifiedBlackbody",
    "NonDissipativePhotosphere",
    "NonDissipativePhotosphere_Deep",
    "GenericFunction",
    "StepFunction",
    "StepFunctionUpper",
    "Sin",
    "DiracDelta",
    "Log_parabola",
    "Exponential_cutoff",
    "PhAbs",
    "TbAbs",
    "WAbs",
    "ZDust",
    "Constant",
    "Line",
    "Quadratic",
    "Cubic",
    "Quartic",
    "Powerlaw",
    "Powerlaw_flux",
    "Powerlaw_Eflux",
    "Cutoff_powerlaw",
    "Cutoff_powerlaw_Ep",
    "Inverse_cutoff_powerlaw",
    "Super_cutoff_powerlaw",
    "SmoothlyBrokenPowerLaw",
    "Broken_powerlaw",
    "Band",
    "Band_grbm",
    "Band_Calderone",
    "DoubleSmoothlyBrokenPowerlaw",
    "Latitude_galactic_diffuse",
    "Gaussian_on_sphere",
    "Asymm_Gaussian_on_sphere",
    "Disk_on_sphere",
    "Ellipse_on_sphere",
    "SpatialTemplate_2D",
    "Power_law_on_sphere",
    "Continuous_injection_diffusion_ellipse",
    "Continuous_injection_diffusion",
    "Continuous_injection_diffusion_legacy",
    "Hermes",
    "Gaussian",
    "Truncated_gaussian",
    "Cauchy",
    "Cosine_Prior",
    "Log_normal",
    "Uniform_prior",
    "Log_uniform_prior",
    "Beta",
    "Gamma",
    "Exponential",
    "Powerlaw_Prior",
    "TemplateModel",
    "TemplateModelFactory",
    "XSPECTableModel",
    "MissingDataFile",
    "DMFitFunction",
    "DMSpectra",
    "functions",
    "list_functions",
]

_functions = {}
for f in _func_names:
    _functions[f] = ("astromodels.functions", f)


_public = {
    # core
    "get_units": ("astromodels.core.units", "get_units"),
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
    "Unpolarized": ("astromodels.core.polarization", "Unpolarized"),
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
    "show_configuration": ("astromodels.utils.configuration", "show_configuration"),
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
