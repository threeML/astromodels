from __future__ import print_function
__author__ = 'giacomov'

# Import IPython display facility, if available. Otherwise,
# create a wrapper which just uses print

try:

    from IPython.display import display

except ImportError:

    def display(*args):
        """
        Mock version of display, used if there is no ipython installed
        """
        print(args)

try:

    from IPython.display import Latex

except ImportError:

    class Latex(object):
        """
        Mock version of the IPython Latex object, used if there is no ipython installed
        """

        def __init__(self, *args, **kwargs):

            pass

        def __repr__(self, *args, **kwargs):

            print("[you need to install IPython to see the Latex representation]")
