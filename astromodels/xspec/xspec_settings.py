from astromodels.xspec import _xspec


def xspec_abund(command_string=None):
    """
    Change/Report the solar abundance table in use. See XSpec manual for help:

    http://heasarc.nasa.gov/xanadu/xspec/manual/XSabund.html

    :param command_string : the command string. If None, returns the current settings, otherwise set to the provided
    settings
    :return: Either none or the current setting
    """

    if command_string is None:

        return _xspec.get_xsabund()

    else:

        _xspec.set_xsabund(command_string)


def xspec_cosmo(command_string=None):
    """
    Define the Cosmology in use within the XSpec models. See Xspec manual for help:

    http://heasarc.nasa.gov/xanadu/xspec/manual/XScosmo.html

    :param command_string : the command string. If None, returns the current settings, otherwise set to the provided
    settings
    :return: Either none or the current setting (H_0, q_0, lambda_0)
    """

    if command_string is None:

        return _xspec.get_xscosmo()

    else:

        _xspec.set_xscosmo(command_string)


def xspec_xsect(command_string=None):
    """
    Change/Report the photoionization cross sections in use for XSpec models. See Xspec manual for help:

    http://heasarc.nasa.gov/xanadu/xspec/manual/XSxsect.html

    :param command_string : the command string. If None, returns the current settings, otherwise set to the provided
    settings
    :return: Either none or the current setting
    """

    if command_string is None:

        return _xspec.get_xsect()

    else:

        _xspec.set_xsect(command_string)
