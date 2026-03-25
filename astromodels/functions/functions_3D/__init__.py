from importlib import import_module
from importlib.util import find_spec

_exports = {
    "Continuous_injection_diffusion": (".functions_3D", []),
    "Continuous_injection_diffusion_ellipse": (".functions_3D", []),
    "Continuous_injection_diffusion_legacy": (".functions_3D", []),
    "GalPropTemplate_3D": (".functions_3D", []),
    "Hermes": (".functions_3D", []),
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
