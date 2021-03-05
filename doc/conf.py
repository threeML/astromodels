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

import sys
import os
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
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The master toctree document.
master_doc = 'index'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '**.ipynb_checkpoints']

source_suffix = ['.rst', '.md']

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#

pygments_style = 'sphinx'

html_theme = 'sphinx_rtd_theme'

html_theme_options = {
    'style_external_links': True,
    # 'vcs_pageview_mode': 'edit',
    'style_nav_header_background': '#0B4BA8',
    #'only_logo': False,
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

version = 'latest'
# The full version, including alpha/beta/rc tags.
release = 'latest'
