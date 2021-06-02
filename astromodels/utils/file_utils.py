# This file contains some defaults, like locations of files, which should not
# change much but benefits anyway of being in one central location

import os
from pathlib import Path

import pkg_resources


def _get_data_file_path(data_file: str) -> Path:
    """
    Returns the absolute path to the required data files.

    :param data_file: relative path to the data file, relative to the astromodels/data path.
    So to get the path to data/dark_matter/gammamc_dif.dat you need to use data_file="dark_matter/gammamc_dif.dat"
    :return: absolute path of the data file
    """

    try:

        file_path = pkg_resources.resource_filename("astromodels", 'data/%s' % data_file)

    except KeyError:

        raise IOError("Could not read or find data file %s. Try reinstalling astromodels. If this does not fix your "
                      "problem, open an issue on github." % (data_file))

    else:

        return Path(file_path).absolute()



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

