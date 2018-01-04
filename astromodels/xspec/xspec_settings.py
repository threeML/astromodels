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


def xspec_cosmo(H0=None,q0=None,lambda_0=None):
    """
    Define the Cosmology in use within the XSpec models. See Xspec manual for help:

    http://heasarc.nasa.gov/xanadu/xspec/manual/XScosmo.html
    
    All parameters can be modified or just a single parameter

    :param H0: the hubble constant
    :param q0:
    :param lambda_0:
    :return: Either none or the current setting (H_0, q_0, lambda_0)
    """

    current_settings = _xspec.get_xscosmo()

    if (H0 is None) and (q0 is None) and (lambda_0 is None):

        return current_settings


    else:

        # ok, we will see what was changed by the used

        user_inputs = [H0, q0, lambda_0]

        for i, current_setting in enumerate(current_settings):

            if user_inputs[i] is None:

                # the user didn't modify this,
                # so lets keep what was already set

                user_inputs[i] = current_setting


        # pass this to xspec

        _xspec.set_xscosmo(*user_inputs)


def xspec_xsect(command_string=None):
    """
    Change/Report the photoionization cross sections in use for XSpec models. See Xspec manual for help:

    http://heasarc.nasa.gov/xanadu/xspec/manual/XSxsect.html

    :param command_string : the command string. If None, returns the current settings, otherwise set to the provided
    settings
    :return: Either none or the current setting
    """

    if command_string is None:

        return _xspec.get_xsxsect()

    else:

        _xspec.set_xsxsect(command_string)
