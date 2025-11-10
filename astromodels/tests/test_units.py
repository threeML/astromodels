from copy import deepcopy
import astropy.units as u
from omegaconf import OmegaConf

from astromodels.core.units import get_units, set_units
from astromodels.utils.configuration import astromodels_config
from astromodels.utils.file_utils import get_path_of_user_config
from astromodels.functions.functions_1D.powerlaws import Powerlaw
from astromodels.sources.point_source import PointSource


def test_config_unit_same():
    for user_config_file in get_path_of_user_config().glob("*.yml"):
        user_conf = OmegaConf.load(user_config_file)
        if "units" in user_conf.keys():
            for k in user_conf["units"].keys():
                assert u.Unit(user_conf["units"][k]) == getattr(get_units(), k)
                assert astromodels_config["units"][k] == user_conf["units"][k]
    for k, v in astromodels_config["units"].items():
        assert u.Unit(v) == getattr(get_units(), k)


def test_set_units():
    mapping = {"energy": u.TeV, "area": u.m**2, "time": u.h, "angle": u.rad}
    ref = deepcopy(get_units())
    not_changed = []
    for k, v in mapping.items():
        if getattr(ref, k) != v:
            set_units(k, v)
        else:
            not_changed.append(k)
    for k, v in mapping.items():
        if k not in not_changed:
            assert getattr(get_units(), k) == v
            assert getattr(ref, k) != getattr(get_units(), k)

    # reset it to the original
    for k, v in mapping.items():
        set_units(k, getattr(ref, k))
