from typing import Optional
import warnings

from omegaconf import OmegaConf
from omegaconf.dictconfig import DictConfig
from rich.tree import Tree

from astromodels.utils import get_path_of_user_config

from .config_structure import Config

# Read the default Config
astromodels_config: Config = OmegaConf.structured(Config)


# now glob the config directory
def update_config_with_user_configs(astromodels_config):
    for user_config_file in get_path_of_user_config().glob("*.yml"):
        _partial_conf = OmegaConf.load(user_config_file)
        if "logging" in _partial_conf.keys():
            if "startup_warnings" in _partial_conf["logging"].keys():
                warnings.warn(
                    "You've provided 'logging.startup_warnings' in "
                    + str(user_config_file)
                    + ". "
                    + "This is deprecated since v2.6.0 - will ignore it"
                )
                del _partial_conf.logging.startup_warnings

        astromodels_config: Config = OmegaConf.merge(astromodels_config, _partial_conf)
    return astromodels_config


astromodels_config = update_config_with_user_configs(astromodels_config)


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
    """Display the current configuration or a sub menu if provided."""

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
