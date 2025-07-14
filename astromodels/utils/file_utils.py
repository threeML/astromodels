# This file contains some defaults, like locations of files, which should not
# change much but benefits anyway of being in one central location

import os
from importlib.resources import files
from pathlib import Path
from typing import Optional, Union

import numpy as np

_custom_config_path = os.environ.get("ASTROMODELS_CONFIG")

copy_if_needed: Optional[bool]

if np.lib.NumpyVersion(np.__version__) >= "2.0.0":
    copy_if_needed = None
else:
    copy_if_needed = False


def _get_data_file_path(data_file: Union[str, Path]) -> Path:
    """Returns the absolute path to the required data files.

    :param data_file: relative path to the data file (str or Path), relative to
    astromodels/data/.
    Example: "dark_matter/gammamc_dif.dat" or Path("dark_matter/gammamc_dif.dat")
    :return: absolute path of the data file as a Path object
    """
    data_file = Path(data_file)

    try:
        resource_path = files("astromodels").joinpath("data", *data_file.parts)

        if not resource_path.is_file():
            raise FileNotFoundError

    except Exception:
        raise IOError(
            f"Could not read or find data file {data_file}. "
            "Try reinstalling astromodels. "
            f"If this does not fix your problem, open an issue on github."
        )

    else:
        return Path(resource_path).resolve()


def get_user_path() -> Path:

    user_path: Path = Path.home() / ".astromodels"

    if not user_path.exists():

        user_path.mkdir(parents=True)

    return user_path


def get_user_data_path() -> Path:

    user_data: Path = get_user_path() / "data"

    # Create it if doesn't exist
    if not user_data.exists():

        user_data.mkdir(parents=True)

    return user_data


def get_path_of_user_config() -> Path:

    if _custom_config_path is not None:

        config_path: Path = Path(_custom_config_path)

    config_path: Path = Path().home() / ".config" / "astromodels"

    if not config_path.exists():

        config_path.mkdir(parents=True)

    return config_path
