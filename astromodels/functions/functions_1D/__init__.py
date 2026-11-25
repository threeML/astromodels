from importlib import import_module
from importlib.util import find_spec

# class name: (file, external dependencies)
_exports = {
    "Blackbody": (".blackbody", []),
    "ModifiedBlackbody": (".blackbody", []),
    "NonDissipativePhotosphere": (".blackbody", []),
    "NonDissipativePhotosphere_Deep": (".blackbody", []),
    "DiracDelta": (".functions", []),
    "Exponential_cutoff": (".functions", []),
    "GenericFunction": (".functions", []),
    "Log_parabola": (".functions", []),
    "Sin": (".functions", []),
    "StepFunction": (".functions", []),
    "StepFunctionUpper": (".functions", []),
    "Synchrotron": (".functions", ["naima"]),
    "Cutoff_powerlaw_flux": (".functions", ["pygsl"]),
    "Band": (".powerlaws", []),
    "Band_Calderone": (".powerlaws", []),
    "Band_grbm": (".powerlaws", []),
    "Broken_powerlaw": (".powerlaws", []),
    "Cutoff_powerlaw": (".powerlaws", []),
    "Cutoff_powerlaw_Ep": (".powerlaws", []),
    "DoubleSmoothlyBrokenPowerlaw": (".powerlaws", []),
    "Inverse_cutoff_powerlaw": (".powerlaws", []),
    "Powerlaw": (".powerlaws", []),
    "Powerlaw_Eflux": (".powerlaws", []),
    "Powerlaw_flux": (".powerlaws", []),
    "SmoothlyBrokenPowerLaw": (".powerlaws", []),
    "Super_cutoff_powerlaw": (".powerlaws", []),
    "APEC": (".apec", ["pyatomdb"]),
    "VAPEC": (".apec", ["pyatomdb"]),
    "Standard_Rv": (".extinction", []),
    "ZDust": (".extinction", []),
    "Constant": (".polynomials", []),
    "Cubic": (".polynomials", []),
    "Line": (".polynomials", []),
    "Quadratic": (".polynomials", []),
    "Quartic": (".polynomials", []),
    "get_polynomial": (".polynomials", []),
    "EBLattenuation": (".absorption", ["ebltable"]),
    "PhAbs": (".absorption", []),
    "WAbs": (".absorption", []),
    "TbAbs": (".absorption", []),
    "has_ebltable": (".absorption", []),
    "has_atomdb": (".apec", []),
    "has_naima": (".functions", []),
    "has_gsl": (".functions", []),
}


def _available(name: str) -> bool:
    _, deps = _exports[name]
    return all(find_spec(d) is not None for d in deps)


# Compute __all__ dynamically based on availability
__all__ = sorted([name for name in _exports if _available(name)])


def __getattr__(name: str):
    # If a name isn’t in __all__, treat as absent
    if name not in __all__:
        # Optional: more specific message
        if name in _exports:
            _, deps = _exports[name]
            missing = [d for d in deps if find_spec(d) is None]
            raise AttributeError(
                f"{name} is unavailable; missing deps: {', '.join(missing)}. "
            )
        raise AttributeError(name)
    mod_path, _ = _exports[name]
    mod = import_module(mod_path, __name__)
    obj = getattr(mod, name)
    globals()[name] = obj
    return obj


def __dir__():
    return sorted(__all__)
