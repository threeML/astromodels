import astropy.units as u
from omegaconf import OmegaConf

from astromodels.core.units import get_units
from astromodels.utils.configuration import astromodels_config
from astromodels.utils.file_utils import get_path_of_user_config


def test_config_unit_same():
    for user_config_file in get_path_of_user_config().glob("*.yml"):
        user_conf = OmegaConf.load(user_config_file)
        if "units" in user_conf.keys():
            for k in user_conf["units"].keys():
                assert u.Unit(user_conf["units"][k]) == getattr(get_units(), k)
                assert astromodels_config["units"][k] == user_conf["units"][k]
    for k, v in astromodels_config["units"].items():
        assert u.Unit(v) == getattr(get_units(), k)
