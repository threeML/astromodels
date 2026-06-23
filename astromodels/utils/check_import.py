import logging
from importlib.util import find_spec

log = logging.getLogger(__name__)


def check_import(module_name, class_name="Function"):
    if find_spec(module_name) is None:
        # only warn here so we can raise custom exceptions or try something else
        log.warning(
            f"{module_name} is required for {class_name} that you are trying to use."
            "Install it and try again."
        )
        return False
    else:
        return True
