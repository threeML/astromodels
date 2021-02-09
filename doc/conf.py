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


DOCS = Path(__file__).parent

# -- Generate API documentation ------------------------------------------------


def run_apidoc(app):
    """Generage API documentation"""
    import better_apidoc

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
            str(DOCS / ".." / "threeML"),
        ]
    )


import astromodels

sys.path.insert(0, os.path.abspath('..'))
sys.path.insert(1, os.path.abspath('../astromodels'))

# Add the path to the C extension
lib_path = os.path.abspath('%s/core' % astromodels.__path__[0])

sys.path.insert(2, lib_path)

# This must work now
import node_ctype


# -- Project information -----------------------------------------------------

project = u'Astromodels'
copyright = u'2016--2020, G.Vianello, J. M. Burgess, N. Di Lalla, N. Omodei, H. Fleischhack'
author = u'G.Vianello'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    'nbsphinx',
    'recommonmark',
    'sphinx.ext.autodoc',
    'sphinx.ext.mathjax',
    'sphinx.ext.viewcode',
#    'sphinx_gallery.gen_gallery',
    'sphinx_gallery.load_style',
    'rtds_actions'
]



napoleon_google_docstring = True
napoleon_use_param = False


# The path where the artifact should be extracted
# Note: this is relative to the conf.py file!
rtds_action_path = "notebooks"
# # The "prefix" used in the `upload-artifact` step of the action
rtds_action_artifact_prefix = "notebooks-for-"


rtds_action_github_repo = "threeML/astromodels"

# # A GitHub personal access token is required, more info below
rtds_action_github_token = os.environ["GITHUB_TOKEN"]


# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The master toctree document.
master_doc = 'index'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '**.ipynb_checkpoints']

source_suffix = ['.rst', '.ipynb']

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#


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

html_theme = 'sphinx_rtd_theme'

# html_theme_options = {
#     'style_external_links': True,
#     # 'vcs_pageview_mode': 'edit',
# #    'style_nav_header_background': '#0B4BA8',
#     # 'only_logo': False,
# }

html_theme_options = {
    'logo_only':False,
    'display_version': False,
    'collapse_navigation': True,
    'navigation_depth': 4,
    'prev_next_buttons_location': 'bottom',  # top and bottom
}


# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

version = 'latest'
# The full version, including alpha/beta/rc tags.
release = 'latest'


def setup(app):
    app.connect("builder-inited", run_apidoc)
