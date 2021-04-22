# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
# import os
# import sys
# sys.path.insert(0, os.path.abspath('.'))


import os
import sys

from pathlib import Path

import mock

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
sys.path.insert(0, os.path.abspath('..'))
#sys.path.insert(0, os.path.abspath('../threeML/classicMLE'))

print(f" current dir {os.getcwd()}")

files = [f for f in os.listdir('.') if os.path.isfile(f)]
for f in files:
    print(f)


DOCS = Path(__file__).parent

# -- Generate API documentation ------------------------------------------------


def run_apidoc(app):
    """Generage API documentation"""
    import better_apidoc
    import pkgutil
    import sys
    import os
    astro_path = os.path.dirname(pkgutil.get_loader("astromodels").get_filename())

    sys.path.insert(0, os.path.abspath('..'))
    sys.path.insert(1, os.path.abspath('../astromodels'))

    # Add the path to the C extension
    #lib_path = os.path.abspath('%s/core' % astromodels.__path__[0])
    lib_path = os.path.abspath('%s/core' % astro_path)

    sys.path.insert(2, lib_path)

    #This must work now
    import node_ctype


    better_apidoc.APP = app
    better_apidoc.main(
        [
            "better-apidoc",
            # "-t",
            # str(docs / "_templates"),
            "--force",
            "--no-toc",
            "--separate",
            "-o",
            str(DOCS / "api"),
            str(DOCS / ".." / "astromodels"),
        ]
    )


# #import astromodels
import pkgutil
astro_path = os.path.dirname(pkgutil.get_loader("astromodels").get_filename())

sys.path.insert(1, os.path.abspath('../astromodels'))

# Add the path to the C extension
#lib_path = os.path.abspath('%s/core' % astromodels.__path__[0])
lib_path = os.path.abspath('%s/core' % astro_path)

sys.path.insert(2, lib_path)




#This must work now
import node_ctype

print(f" current dir {os.getcwd()}")
files = [f for f in os.listdir('.') if os.path.isfile(f)]
for f in files:
    print(f)
    
# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'rtds_action',
    'nbsphinx',
    'recommonmark',
    'sphinx.ext.autodoc',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
    'sphinx.ext.githubpages',
    'sphinx.ext.napoleon',
    'sphinx_gallery.load_style',
    
]



napoleon_google_docstring = True
napoleon_use_param = False


# The path where the artifact should be extracted
# Note: this is relative to the conf.py file!

rtds_action_path = "notebooks"

# # The "prefix" used in the `upload-artifact` step of the action
rtds_action_artifact_prefix = "notebooks-for-"


rtds_action_github_repo = "astromodels/astromodels"

# # A GitHub personal access token is required, more info below
rtds_action_github_token = os.environ["GITHUB_TOKEN"]

rtds_action_error_if_missing = True


# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

source_suffix = ['.rst', '.ipynb']

# The master toctree document.
master_doc = 'index'



# -- Project information -----------------------------------------------------

project = u'Astromodels'
copyright = u'2016--2021, G.Vianello, J. M. Burgess, N. Di Lalla, N. Omodei, H. Fleischhack'
author = u'G.Vianello'




# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
language = None

pygments_style = 'none'

# Create local pygments copies
# Previously used: https://github.com/richleland/pygments-css
# But do not want to depend on some random repository
from pygments.formatters import HtmlFormatter  # noqa: E402
from pygments.styles import get_all_styles  # noqa: E402
path = os.path.join('_static', 'pygments')
if not os.path.isdir(path):
    os.mkdir(path)
for style in get_all_styles():
    path = os.path.join('_static', 'pygments', style + '.css')
    if os.path.isfile(path):
        continue
    with open(path, 'w') as f:
        f.write(HtmlFormatter(style=style).get_style_defs('.highlight'))




# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', '**.ipynb_checkpoints', 'md/*.md']

pygments_style = 'sphinx'

html_theme = 'sphinx_rtd_theme'


html_theme_options = {
    'logo_only':False,
    'display_version': False,
    'collapse_navigation': True,
    'navigation_depth': 4,
    'prev_next_buttons_location': 'bottom',  # top and bottom
}


html_logo = "media/transp_logo.png"
html_show_sourcelink = False
html_favicon = "media/favicon.ico"

autosectionlabel_prefix_document = True

html_static_path = ['_static']

version = 'latest'
# The full version, including alpha/beta/rc tags.
release = 'latest'

print("Done.")

# def setup(app):
#     app.connect("builder-inited", run_apidoc)
