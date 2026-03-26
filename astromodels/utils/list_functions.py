from __future__ import annotations

import collections
import inspect
import sys
from importlib import import_module
from importlib.util import find_spec
from typing import Any, Dict, List, Optional

from astromodels.utils.exceptions import UnknownFunction
from astromodels.utils.table import dict_to_table

dict_of_available_functions = None


def _iter_registry() -> Dict[str, Any]:
    """
    Return a shallow copy of the registry: public name -> class.
    """
    from astromodels.functions.function import iter_registered_functions  # type: ignore

    try:
        return iter_registered_functions()
    except Exception:
        # Fallback if helper is not present; try the internal name
        from astromodels.functions.function import _FUNCTION_REGISTRY  # type: ignore

        return dict(_FUNCTION_REGISTRY)


def _description_for(cls: Any) -> str:
    """
    Extract a human-friendly description for a Function subclass.
    """
    # Prefer the parsed YAML definition if present
    try:
        desc = getattr(cls, "_function_definition", {}).get("description", "")
        if desc:
            return desc
    except Exception:
        pass

    # Fallback to first line of the docstring
    doc = getattr(cls, "__doc__", "") or ""
    if doc.strip():
        return doc.strip().splitlines()[0]
    return ""


def list_functions(return_dict: bool = False):
    """
    Returns all registered functions. Even unavailable ones

    :param return_dict: return it as a dict instead of an astropy.table.Table.
        Defaults to false.
    :type return_dict: bool

    :returns: astropy.table.Table or dict

    """
    functions_and_descriptions: Dict[str, Dict[str, str]] = {}

    reg = _iter_registry()
    for name, cls in reg.items():
        functions_and_descriptions[name] = {"Description": _description_for(cls)}

    public_funcs: List[str] = []
    try:
        import astromodels.functions as fpkg  # noqa: WPS433

        public_funcs = list(getattr(fpkg, "__all__", []))
    except Exception:
        public_funcs = []

    public_xs: List[str] = []
    xs_mod = sys.modules.get("astromodels.xspec.factory")
    if xs_mod is not None:
        try:
            public_xs = list(getattr(xs_mod, "__all__", []))
        except Exception:
            public_xs = []

    registered_names = set(reg.keys())
    skip_names = {
        "list_functions",
        "has_atomdb",
        "FunctionMeta",
        "Function1D",
        "Function2D",
        "Function3D",
        "get_polynomial",
    }

    for name in public_funcs + public_xs:
        if name in registered_names or name in skip_names:
            continue
        functions_and_descriptions.setdefault(
            name, {"Description": "(unavailable: not imported)"}
        )

    ordered = collections.OrderedDict(sorted(functions_and_descriptions.items()))
    if return_dict:
        return ordered
    return dict_to_table(ordered)


def get_function_class(name: str) -> Any:
    """
    Resolve a function class by name with minimal side effects.

    Resolution order:
    1) Registry (covers built-ins, XS_* created by xspec.factory, runtime classes).
    2) If name does not start with 'XS_', trigger lazy import of that symbol from
       astromodels.functions, then re-check the registry.
    3) If name starts with 'XS_', look it up only if astromodels.xspec.factory is
       already imported (no import here); otherwise fail.

    :param name: name of the function
    :type name: str

    :returns: astromodels.functions.function.Function
    """
    reg = _iter_registry()
    cls = reg.get(name)
    if cls is not None:
        return cls

    if not name.startswith("XS_"):
        try:
            fpkg = import_module("astromodels.functions")
            getattr(fpkg, name)  # triggers lazy import/definition
        except Exception as e:
            raise UnknownFunction from e

        reg = _iter_registry()
        cls = reg.get(name)
        if cls is not None:
            return cls

        try:
            return getattr(fpkg, name)
        except Exception as e:
            raise UnknownFunction from e
    xs_mod: Optional[Any] = sys.modules.get("astromodels.xspec.factory")
    if xs_mod is not None:
        try:
            return getattr(xs_mod, name)
        except AttributeError:
            # Might be re-exported at astromodels.xspec as a convenience
            xpkg = sys.modules.get("astromodels.xspec")
            if xpkg is not None:
                try:
                    return getattr(xpkg, name)
                except AttributeError:
                    pass

    # Not found
    raise UnknownFunction


def list_available_functions(return_dict: bool = False):
    """
    Get a table of all available Functions


    :param return_dict: return it as a dict instead of an astropy.table.Table.
        Defaults to false.
    :type return_dict: bool

    :returns: astropy.table.Table or dict
    """
    global dict_of_available_functions
    if dict_of_available_functions is None:
        import astromodels.functions as fpkg
        from astromodels.functions.function import Function

        skip_names: set[str] = {
            "list_functions",
            "has_atomdb",
            "FunctionMeta",
            "Function1D",
            "Function2D",
            "Function3D",
            "get_polynomial",
        }

        functions_and_descriptions: Dict[str, Dict[str, str]] = {}

        public: List[str] = list(getattr(fpkg, "__all__", []))
        if find_spec("xspec") is not None:
            import astromodels.xspec.factory as xspec_fact

            public.extend(list(getattr(xspec_fact, "__all__", [])))

        for name in public:
            if name in skip_names:
                continue

            try:
                if not name.startswith("XS_"):
                    obj = getattr(fpkg, name)
                else:
                    obj = getattr(xspec_fact, name)

            except Exception:
                continue

            if inspect.isclass(obj):
                if issubclass(obj, Function) and obj is not Function:
                    desc = None
                    try:
                        # Many astromodels function classes carry this dict
                        desc = getattr(obj, "_function_definition", {}).get(
                            "description"
                        )
                    except Exception:
                        desc = None
                    functions_and_descriptions[name] = {"Description": desc}

        ordered = dict(sorted(functions_and_descriptions.items()))
        dict_of_available_functions = ordered.copy()
    else:
        ordered = dict_of_available_functions

    if not return_dict:
        return dict_to_table(ordered)
    return ordered


def list_function_names() -> List[str]:
    """
    Return the list of registered function names (sorted).
    """
    reg = _iter_registry()
    return sorted(reg.keys())


def list_available_function_names() -> List[str]:
    """
    Get only the names without a description from the available Functions
    """
    if dict_of_available_functions is None:
        _ = list_available_functions(return_dict=True)
    return list(dict_of_available_functions.keys())
