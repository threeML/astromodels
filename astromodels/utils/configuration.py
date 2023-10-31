from pathlib import Path
from typing import Any, Dict, Optional

from omegaconf import OmegaConf
from omegaconf.dictconfig import DictConfig
from rich.tree import Tree

from astromodels.utils import get_path_of_user_config

from .config_structure import Config

# Read the default Config
astromodels_config: Config = OmegaConf.structured(Config)

# now glob the config directory

for user_config_file in get_path_of_user_config().glob("*.yml"):

    _partial_conf = OmegaConf.load(user_config_file)

    astromodels_config: Config = OmegaConf.merge(astromodels_config, _partial_conf)


def recurse_dict(d, tree):

    for k, v in d.items():

        if (type(v) is dict) or isinstance(v, DictConfig):

            branch = tree.add(
                k, guide_style="bold medium_orchid", style="bold medium_orchid"
            )

            recurse_dict(v, branch)

        else:

            tree.add(
                f"{k}: [blink cornflower_blue]{v}",
                guide_style="medium_spring_green",
                style="medium_spring_green",
            )

    return


def show_configuration(sub_menu: Optional[str] = None):
    """
    display the current configuration or a sub menu if
    provided
    """

    if sub_menu is None:

        tree = Tree(
            "config",
            guide_style="bold medium_orchid",
            style="bold medium_orchid",
        )

        recurse_dict(astromodels_config, tree)

    else:

        if sub_menu in astromodels_config:

            tree = Tree(
                "config",
                guide_style="bold medium_orchid",
                style="bold medium_orchid",
            )

            recurse_dict(astromodels_config[sub_menu], tree)

        else:

            msg = f"{sub_menu} is not in the astromodels configuration"

            raise AssertionError(msg)

    return tree
