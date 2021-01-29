import os
import re
import sys

# adapted from thepipe
# https://github.com/tamasgal/thepipe/


ATTRIBUTES = dict(
    list(
        zip([
            'bold', 'dark', '', 'underline', 'blink', '', 'reverse',
            'concealed'
        ], list(range(1, 9)))))
del ATTRIBUTES['']

ATTRIBUTES_RE = r'\033\[(?:%s)m' % '|'  \
                .join(['%d' % v for v in ATTRIBUTES.values()])

HIGHLIGHTS = dict(
    list(
        zip([
            'on_grey', 'on_red', 'on_green', 'on_yellow', 'on_blue',
            'on_magenta', 'on_cyan', 'on_white'
        ], list(range(40, 48)))))

HIGHLIGHTS_RE = r'\033\[(?:%s)m' % '|'  \
                .join(['%d' % v for v in HIGHLIGHTS.values()])

COLORS = dict(
    list(
        zip([
            'grey',
            'red',
            'green',
            'yellow',
            'blue',
            'magenta',
            'cyan',
            'white',
        ], list(range(30, 38)))))

COLORS_RE = r'\033\[(?:%s)m' % '|'.join(['%d' % v for v in COLORS.values()])

RESET = r'\033[0m'
RESET_RE = r'\033\[0m'


def colored(text, color=None, on_color=None, attrs=None, ansi_code=None):
    """Colorize text, while stripping nested ANSI color sequences.
    Author:  Konstantin Lepa <konstantin.lepa@gmail.com> / termcolor
    Available text colors:
        red, green, yellow, blue, magenta, cyan, white.
    Available text highlights:
        on_red, on_green, on_yellow, on_blue, on_magenta, on_cyan, on_white.
    Available attributes:
        bold, dark, underline, blink, reverse, concealed.
    Example:
        colored('Hello, World!', 'red', 'on_grey', ['blue', 'blink'])
        colored('Hello, World!', 'green')
    """
    if os.getenv('ANSI_COLORS_DISABLED') is None:
        if ansi_code is not None:
            return "\033[38;5;{}m{}\033[0m".format(ansi_code, text)
        fmt_str = '\033[%dm%s'
        if color is not None:
            text = re.sub(COLORS_RE + '(.*?)' + RESET_RE, r'\1', text)
            text = fmt_str % (COLORS[color], text)
        if on_color is not None:
            text = re.sub(HIGHLIGHTS_RE + '(.*?)' + RESET_RE, r'\1', text)
            text = fmt_str % (HIGHLIGHTS[on_color], text)
        if attrs is not None:
            text = re.sub(ATTRIBUTES_RE + '(.*?)' + RESET_RE, r'\1', text)
            for attr in attrs:
                text = fmt_str % (ATTRIBUTES[attr], text)
        return text + RESET
    else:
        return text


def cprint(text, color=None, on_color=None, attrs=None):
    """Print colorize text.
    Author:  Konstantin Lepa <konstantin.lepa@gmail.com> / termcolor
    It accepts arguments of print function.
    """
    print((colored(text, color, on_color, attrs)))


def isnotebook():
    """Check if running within a Jupyter notebook"""
    try:
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell':
            return True  # Jupyter notebook or qtconsole
        elif shell == 'TerminalInteractiveShell':
            return False  # Terminal running IPython
        else:
            return False  # Other type (?)
    except NameError:
        return False


def supports_color():
    """Checks if the terminal supports color."""
    if isnotebook():
        return True
    supported_platform = sys.platform != 'win32' or 'ANSICON' in os.environ
    is_a_tty = hasattr(sys.stdout, 'isatty') and sys.stdout.isatty()

    if not supported_platform or not is_a_tty:
        return False

    return True
