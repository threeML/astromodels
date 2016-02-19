__author__ = 'giacomov'

import astropy.table


def dict_to_table(dictionary):

    assert len(dictionary.values()) > 0,"Dictionary cannot be empty"

    # Create an empty table

    table = Table()

    # Add the names as first column

    table['name'] = dictionary.keys()

    # Now add all other properties

    # Use the first parameter as prototype

    prototype = dictionary.values()[0]

    column_names = prototype.keys()

    for column_name in column_names:

        table[column_name] = map(lambda x: x[column_name], dictionary.values())

    return table


# A hack on the astropy Table class to make its output
# more appealing, especially when in the Ipython notebook

class Table(astropy.table.Table):
    """
    Wrapper around the astropy table to remove some useless clutter (like the format of each column)
    """

    def _base_repr_(self, html=False, show_name=True, **kwargs):
        """
        Override the method in the astropy.Table class
        to avoid displaying the description, and the format
        of the columns
        """

        table_id = 'table{id}'.format(id=id(self))

        data_lines, outs = self.formatter._pformat_table(self,
                                                         tableid=table_id, html=html, max_width=(-1 if html else None),
                                                         show_name=show_name, show_unit=None, show_dtype=False)

        out = '\n'.join(data_lines)

        # if astropy.table.six.PY2 and isinstance(out, astropy.table.six.text_type):
        #    out = out.encode('utf-8')

        return out


class NumericMatrix(Table):

    def _base_repr_(self, html=False):
        return super(NumericMatrix, self)._base_repr_(html, False)

