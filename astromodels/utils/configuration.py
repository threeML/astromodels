# This file contains some defaults, like locations of files, which should not
# change much but benefits anyway of being in one central location

from pathlib import Path


def get_user_path():

    user_path = Path.home() / ".astromodels"

    if not user_path.exists():

        user_path.mkdir(parents=True)

    return user_path
        
    
def get_user_data_path():

    user_data = get_user_path() / "data"

    # Create it if doesn't exist
    if not user_data.exists():

        user_data.mkdir(parents=True)


    return user_data

