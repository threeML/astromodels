import pandas as pd


def long_path_formatter(line, max_width=pd.get_option('max_colwidth')):
    """
    If a path is longer than max_width, it substitute it with the first and last element,
    joined by "...". For example 'this.is.a.long.path.which.we.want.to.shorten' becomes
    'this...shorten'

    :param line:
    :param max_width:
    :return:
    """

    if len(line) > max_width:

        tokens = line.split(".")
        trial1 = "%s...%s" % (tokens[0], tokens[-1])

        if len(trial1) > max_width:

            return "...%s" %(tokens[-1][-1:-(max_width-3)])

        else:

            return trial1

    else:

        return line
