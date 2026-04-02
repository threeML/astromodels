__author__ = "giacomov"
from importlib.util import find_spec

# Import IPython display facility, if available. Otherwise,
# create a wrapper which just uses print
if find_spec("IPython") is not None:
    from IPython.display import Latex, display
else:

    def display(*args):
        """Mock version of display, used if there is no ipython installed."""
        print(args)

    class Latex:
        """Mock version of the IPython Latex object, used if there is no
        ipython installed."""

        def __init__(self, *args, **kwargs):

            pass

        def __repr__(self, *args, **kwargs):

            print("[you need to install IPython to see the Latex representation]")
