from importlib.util import find_spec


def check_import(module_name, class_name="Function"):
    if find_spec(module_name) is None:
        raise RuntimeError(
            f"{module_name} is required for {class_name} that you are trying to use."
            "Install it and try again."
        )
    else:
        return True
