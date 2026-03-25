import collections
import inspect
from typing import Dict, List


def list_functions(include_unavailable: bool = False):
    """
    Gather all public function classes exposed by `astromodels.functions`,
    collect their descriptions, and return a formatted table (as before).

    Parameters
    ----------
    include_unavailable : bool
        If True, include names in the public API that cannot be imported due
        to missing optional dependencies, with a placeholder description.

    Returns
    -------
    Any
        The same table structure produced by `dict_to_table(OrderedDict)`
        previously. If `dict_to_table` is not available, returns the
        `OrderedDict` directly as a fallback.
    """
    # Lazy import to avoid any early import cycles
    import astromodels.functions as fpkg

    # Import the Function base class from its canonical location
    from astromodels.functions.function import Function
    from astromodels.utils.table import dict_to_table

    # Names we should not try to introspect to avoid recursion / non-class entries
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
    for name in public:
        if name in skip_names:
            continue

        try:
            obj = getattr(fpkg, name)
        except Exception:
            if include_unavailable:
                functions_and_descriptions[name] = {
                    "Description": "(unavailable: missing optional dependency)"
                }
            continue

        if inspect.isclass(obj):
            if issubclass(obj, Function) and obj is not Function:
                desc = None
                try:
                    # Many astromodels function classes carry this dict
                    desc = getattr(obj, "_function_definition", {}).get("description")
                except Exception:
                    desc = None
                functions_and_descriptions[name] = {"Description": desc}

    ordered = collections.OrderedDict(sorted(functions_and_descriptions.items()))
    return dict_to_table(ordered)
