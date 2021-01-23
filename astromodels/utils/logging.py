

# If threeml is not installed, we create our own log,
# otherwise, just print to the 3ML one!


import logging
import logging.handlers as handlers
import sys
from pathlib import Path
from typing import Dict, Optional

import colorama
from colorama import Back, Fore, Style

try:

    from threeML.config.config import threeML_config

    has_threeml = True

except:

    has_threeml = False


colorama.deinit()
colorama.init(strip=False)
# set up the console logging


def get_path_of_log_dir():

    # we use the 3ML log path to simplify things
    # a more clever solution could be found

    if has_threeml:

        user_log: Path = Path(threeML_config["logging"]["path"])

    else:

        user_log = Path().home() / ".astromodels" / "log"

    # Create it if doesn't exist
    if not user_log.exists():

        user_log.mkdir(parents=True)

    return user_log


_log_file_names = ["usr.log", "dev.log"]


def get_path_of_log_file(log_file: str) -> Path:
    """
    returns the path of the log files
    """
    assert log_file in _log_file_names, f"{log_file} is not one of {_log_file_names}"

    return get_path_of_log_dir() / log_file


class ColoredFormatter(logging.Formatter):
    """
    Colored log formatter.
    """

    def __init__(
        self, *args, colors: Optional[Dict[str, str]] = None, **kwargs
    ) -> None:
        """Initialize the formatter with specified format strings."""

        super().__init__(*args, **kwargs)

        self.colors = colors if colors else {}

    def format(self, record) -> str:
        """Format the specified record as text."""

        record.color = self.colors.get(record.levelname, "")
        record.reset = Style.RESET_ALL

        return super().format(record)


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

# now set up the console logger

_console_formatter = ColoredFormatter(
    "{color} {levelname:8} {reset}| {color} {message} {reset}",
    style="{",
    datefmt="%H:%M:%S",
    colors={
        "DEBUG": Fore.CYAN,
        "INFO": Fore.GREEN + Style.BRIGHT,
        "WARNING": Fore.YELLOW + Style.DIM,
        "ERROR": Fore.RED + Style.BRIGHT,
        "CRITICAL": Fore.RED + Back.WHITE + Style.BRIGHT,
    },
)

astromodels_console_log_handler = logging.StreamHandler(sys.stdout)
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
    log.setLevel(logging.DEBUG)

    # add the handlers

    log.addHandler(astromodels_dev_log_handler)

    log.addHandler(astromodels_console_log_handler)

    log.addHandler(astromodels_usr_log_handler)

    # we do not want to duplicate teh messages in the parents
    log.propagate = False

    return log
