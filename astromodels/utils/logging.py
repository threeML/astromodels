# If threeml is not installed, we create our own log,
# otherwise, just print to the 3ML one!


import logging
import logging.handlers as handlers
from pathlib import Path

from astromodels.utils.configuration import astromodels_config
from rich.console import Console
from rich.logging import RichHandler
from rich.theme import Theme

from .file_utils import _get_data_file_path

try:

    from threeML.config.config import threeML_config

    has_threeml = True

except ImportError:

    has_threeml = False


DEBUG_NODE_LEVEL = 9
logging.addLevelName(DEBUG_NODE_LEVEL, "DEBUG_NODE")


def debug_node(self, message, *args, **kws):
    if self.isEnabledFor(DEBUG_NODE_LEVEL):
        # Yes, logger takes its '*args' as 'args'.
        self._log(DEBUG_NODE_LEVEL, message, args, **kws)


logging.Logger.debug_node = debug_node


def get_path_of_log_dir():

    # we use the 3ML log path to simplify things
    # a more clever solution could be found

    if has_threeml:

        user_log: Path = Path(threeML_config.logging.path).expanduser()

    else:

        user_log: Path = Path(astromodels_config.logging.path).expanduser()

    # Create it if doesn't exist
    if not user_log.exists():

        user_log.mkdir(parents=True)

    return user_log


_log_file_names = ["usr.log", "dev.log"]


def get_path_of_log_file(log_file: str) -> Path:
    """
    returns the path of the log files
    """
    assert (
        log_file in _log_file_names
    ), f"{log_file} is not one of {_log_file_names}"

    return get_path_of_log_dir() / log_file


class LogFilter(object):
    def __init__(self, level):
        self.__level = level

    def filter(self, logRecord):
        return logRecord.levelno != self.__level


# now create the developer handler that rotates every day and keeps
# 10 days worth of backup
astromodels_dev_log_handler = handlers.TimedRotatingFileHandler(
    get_path_of_log_file("dev.log"), when="D", interval=1, backupCount=10
)

# lots of info written out

_dev_formatter = logging.Formatter(
    "%(asctime)s | %(name)s | %(levelname)s| %(funcName)s | %(lineno)d | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

astromodels_dev_log_handler.setFormatter(_dev_formatter)
astromodels_dev_log_handler.setLevel(logging.DEBUG)
# now set up the usr log which will save the info

astromodels_usr_log_handler = handlers.TimedRotatingFileHandler(
    get_path_of_log_file("usr.log"), when="D", interval=1, backupCount=10
)

astromodels_usr_log_handler.setLevel(logging.INFO)

# lots of info written out
_usr_formatter = logging.Formatter(
    "%(asctime)s | %(levelname)s | %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
)

astromodels_usr_log_handler.setFormatter(_usr_formatter)


_console_formatter = logging.Formatter(
    " %(message)s",
    datefmt="%H:%M:%S",
)


_theme = {}

# Banner
_theme["h1"] = "deep_sky_blue3"
_theme["status.spinner"] = "cyan2"
_theme["status.text"] = "deep_sky_blue4"
_theme["repr.filename"] = "blue"
_theme["repr.number"] = "white"
_theme["repr.path"] = "grey37"
_theme["repr.str"] = "grey37"
_theme["repr.tag_name"] = "white"
_theme["repr.url"] = "not bold not italic underline grey84"
_theme["log.time"] = "green1"
_theme["log.message"] = f"{astromodels_config.logging.message_style}"
_theme["logging.level.debug"] = f"{astromodels_config.logging.debug_style}"
_theme["logging.level.error"] = f"{astromodels_config.logging.error_style}"
_theme["logging.level.info"] = f"{astromodels_config.logging.info_style}"
_theme["logging.level.warning"] = f"{astromodels_config.logging.warn_style}"
_theme["logging.level.degub_node"] = "light_goldenrod1"


# mytheme = Theme().read(_get_data_file_path("log_theme.ini"))
mytheme = Theme(_theme)
console = Console(theme=mytheme)


astromodels_console_log_handler = RichHandler(
    level="INFO", rich_tracebacks=True, markup=True, console=console
)
astromodels_console_log_handler.setFormatter(_console_formatter)
astromodels_console_log_handler.setLevel("INFO")

warning_filter = LogFilter(logging.WARNING)


def silence_warnings():
    """
    supress warning messages in console and file usr logs
    """

    astromodels_usr_log_handler.addFilter(warning_filter)
    astromodels_console_log_handler.addFilter(warning_filter)


def activate_warnings():
    """
    supress warning messages in console and file usr logs
    """

    astromodels_usr_log_handler.removeFilter(warning_filter)
    astromodels_console_log_handler.removeFilter(warning_filter)


def update_logging_level(level):

    astromodels_console_log_handler.setLevel(level)


def setup_logger(name):

    # A logger with name name will be created
    # and then add it to the print stream
    log = logging.getLogger(name)

    # this must be set to allow debug messages through
    log.setLevel(DEBUG_NODE_LEVEL)

    # add the handlers

    log.addHandler(astromodels_dev_log_handler)

    log.addHandler(astromodels_console_log_handler)

    log.addHandler(astromodels_usr_log_handler)

    # we do not want to duplicate teh messages in the parents
    log.propagate = False

    return log
