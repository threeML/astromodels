__author__ = "giacomov"

from astropy.table import Table


def dict_to_table(dictionary, list_of_keys=None):
    """Return a table representing the dictionary.

    :param dictionary: the dictionary to represent
    :param list_of_keys: optionally, only the keys in this list will be
        inserted in the table
    :return: a Table instance
    """

    # assert len(dictionary.values()) > 0, "Dictionary cannot be empty"

    # Create an empty table

    table = Table()

    # If the dictionary is not empty, fill the table

    if len(dictionary) > 0:

        # Add the names as first column

        table["name"] = list(dictionary.keys())

        # Now add all other properties

        # Use the first parameter as prototype

        prototype = list(dictionary.values())[0]

        column_names = list(prototype.keys())

        # If we have a white list for the columns, use it

        if list_of_keys is not None:

            column_names = [key for key in column_names if key in list_of_keys]

        # Fill the table

        for column_name in column_names:

            table[column_name] = [x[column_name] for x in list(dictionary.values())]

    return table
