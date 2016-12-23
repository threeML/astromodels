import pkg_resources
import os


def _get_data_file_path(data_file):
    """
    Returns the absolute path to the required data files.

    :param data_file: relative path to the data file, relative to the astromodels/data path.
    So to get the path to data/dark_matter/gammamc_dif.dat you need to use data_file="dark_matter/gammamc_dif.dat"
    :return: absolute path of the data file
    """

    try:

        file_path = pkg_resources.resource_filename("astromodels", 'data/%s' % data_file)

    except KeyError:

        raise IOError("Could not read or find data file %s. Try reinstalling astromodels. If this does not fix your "
                      "problem, open an issue on github." % (data_file))

    else:

        return os.path.abspath(file_path)
