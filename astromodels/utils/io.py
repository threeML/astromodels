__author__ = 'giacomov'

# Import IPython display facility, if available. Otherwise,
# create a wrapper which just uses print

try:

    from IPython.display import display

except ImportError:

    def display(*args):

        print(args)

