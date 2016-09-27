from ast import parse


def is_valid_variable_name(string_to_check):
    """
    Returns whether the provided name is a valid variable name in Python

    :param string_to_check: the string to be checked
    :return: True or False
    """

    try:

        parse('{} = None'.format(string_to_check))
        return True

    except (SyntaxError, ValueError, TypeError):

        return False
