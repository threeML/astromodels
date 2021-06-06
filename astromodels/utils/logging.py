

# If threeml is not installed, we create our own log,
# otherwise, just print to the 3ML one!


import logging
import logging.handlers as handlers
import sys
from hashlib import sha256
from pathlib import Path

from astromodels.utils.configuration import astromodels_config

from .sys_tools import colored, supports_color

try:

    from threeML.config.config import threeML_config

    has_threeml = True

except:

    has_threeml = False


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
    assert log_file in _log_file_names, f"{log_file} is not one of {_log_file_names}"

    return get_path_of_log_dir() / log_file


DEFAULT_LOG_COLORS = {
    'DEBUG': 'blue',
    'INFO': 'green',
    'WARNING': 'purple',
    'ERROR': 'red',
    'CRITICAL': 'bold_red',

}


def esc(*x):
    """Create escaped code from format code"""
    return '\033[' + ';'.join(x) + 'm'


# The following coloured log logic is from
# https://github.com/borntyping/python-colorlog
# I dropped some features and removed the Python 2.7 compatibility

ESCAPE_CODES = {'reset': esc('0'), 'bold': esc('01'), 'thin': esc('02')}

COLORS = ['black', 'red', 'green', 'yellow', 'blue', 'purple', 'cyan', 'white']

PREFIXES = [
    # Foreground without prefix
    ('3', ''),
    ('01;3', 'bold_'),
    ('02;3', 'thin_'),

    # Foreground with fg_ prefix
    ('3', 'fg_'),
    ('01;3', 'fg_bold_'),
    ('02;3', 'fg_thin_'),

    # Background with bg_ prefix - bold/light works differently
    ('4', 'bg_'),
    ('10', 'bg_bold_'),
]

for _prefix, _prefix_name in PREFIXES:
    for _code, _name in enumerate(COLORS):
        ESCAPE_CODES[_prefix_name + _name] = esc(_prefix + str(_code))


def parse_colors(sequence):
    """Return escape codes from a color sequence."""
    return ''.join(ESCAPE_CODES[n] for n in sequence.split(',') if n)


class ColoredRecord:
    """
    Wraps a LogRecord, adding named escape codes to the internal dict.
    The internal dict is used when formatting the message (by the PercentStyle,
    StrFormatStyle, and StringTemplateStyle classes).
    """

    def __init__(self, record):
        """Add attributes from the escape_codes dict and the record."""
        self.__dict__.update(ESCAPE_CODES)
        self.__dict__.update(record.__dict__)
        self.__record = record

    def __getattr__(self, name):
        return getattr(self.__record, name)


class ColoredFormatter(logging.Formatter):
    """
    A formatter that allows colors to be placed in the format string.
    Intended to help in creating more readable logging output.
    Based on https://github.com/borntyping/python-colorlog
    """

    def __init__(self,
                 fmt,
                 datefmt=None,
                 style='%',
                 log_colors=None,
                 reset=True,
                 secondary_log_colors=None):
        """
        Set the format and colors the ColouredFormatter will use.
        The ``fmt``, ``datefmt`` and ``style`` args are passed on to the
        ``logging.Formatter`` constructor.
        The ``secondary_log_colors`` argument can be used to create additional
        ``log_color`` attributes. Each key in the dictionary will set
        ``{key}_log_color``, using the value to select from a different
        ``log_colors`` set.
        :Parameters:
        - fmt (str): The format string to use
        - datefmt (str): A format string for the date
        - log_colors (dict):
            A mapping of log level names to color names
        - reset (bool):
            Implicitly append a color reset to all records unless False
        - style ('%' or '{' or '$'):
            The format style to use. (*No meaning prior to Python 3.2.*)
        - secondary_log_colors (dict):
            Map secondary ``log_color`` attributes. (*New in version 2.6.*)
        """
        super(ColoredFormatter, self).__init__(fmt, datefmt, style)

        self.log_colors = (log_colors
                           if log_colors is not None else DEFAULT_LOG_COLORS)
        self.secondary_log_colors = secondary_log_colors
        self.reset = reset

    def format(self, record):
        """Format a message from a record object."""
        record = ColoredRecord(record)
        record.log_color = escape_codes(self.log_colors, record.levelname)

        if self.secondary_log_colors:
            for name, log_colors in self.secondary_log_colors.items():
                color = escape_codes(log_colors, record.levelname)
                setattr(record, name + '_log_color', color)

        message = super(ColoredFormatter, self).format(record)

        if self.reset and not message.endswith(ESCAPE_CODES['reset']):
            message += ESCAPE_CODES['reset']

        return message


def escape_codes(log_colors, level_name):
    """Return escape codes from a ``log_colors`` dict."""
    return parse_colors(log_colors.get(level_name, ""))


def hash_coloured(text):
    """Return a ANSI coloured text based on its hash"""
    ansi_code = int(sha256(text.encode('utf-8')).hexdigest(), 16) % 230
    return colored(text, ansi_code=ansi_code)


def hash_coloured_escapes(text):
    """Return the ANSI hash colour prefix and suffix for a given text"""
    ansi_code = int(sha256(text.encode('utf-8')).hexdigest(), 16) % 230
    prefix, suffix = colored('SPLIT', ansi_code=ansi_code).split('SPLIT')
    return prefix, suffix



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
name = "test"
if supports_color():
    prefix_1, suffix = hash_coloured_escapes(name)
    prefix_2, _ = hash_coloured_escapes(name + 'salt')
else:
    prefix_1, prefix_2, suffix = ('', '', '')

date_str = ''

_console_formatter = ColoredFormatter('[%(log_color)s%(levelname)-8s%(reset)s]'
                                      '%(log_color)s %(message)s',
                                      datefmt="%H:%M:%S",

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
