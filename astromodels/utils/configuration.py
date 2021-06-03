from omegaconf import OmegaConf

from astromodels.utils import get_path_of_user_config

from .config_structure import Config

# Read the default Config
astromodels_config: Config = OmegaConf.structured(Config)

# now glob the config directory

for user_config_file in get_path_of_user_config().glob("*.yml"):

    _partial_conf = OmegaConf.load(user_config_file)

    astromodels_config: Config = OmegaConf.merge(astromodels_config, _partial_conf)
