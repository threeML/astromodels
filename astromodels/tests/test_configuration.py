import os

import yaml
from omegaconf import OmegaConf

from astromodels.utils.configuration import update_config_with_user_configs
from astromodels.utils.config_structure import Config


def test_user_configuration(tmp_path):
    original_config_path = os.environ.get("ASTROMODELS_CONFIG")
    os.environ["ASTROMODELS_CONFIG"] = str(tmp_path)
    try:

        dummy_config = OmegaConf.structured(Config)

        configs = [
            {"logging": {"usr": "off", "startup_warnings": "off"}},
        ]

        for i, c in enumerate(configs):
            path = tmp_path / f"conf_{i}.yml"

            with path.open("w") as f:
                yaml.dump(stream=f, data=c, Dumper=yaml.SafeDumper)

        dummy_config = update_config_with_user_configs(dummy_config)
    except Exception as e:
        raise e
    finally:
        original_config_path = (
            "" if original_config_path is None else original_config_path
        )
        os.environ["ASTROMODELS_CONFIG"] = original_config_path
