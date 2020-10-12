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
import os
import sys
import sphinx_rtd_theme
sys.path.insert(0, os.path.abspath('../'))


# -- Project information -----------------------------------------------------

project = 'ConvTrack'
copyright = '2020, James Dean'
author = 'James Dean'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'docs.ext.*') or your custom
# ones.
extensions = [
    'docs.ext.autodoc',
    'docs.ext.doctest',
    'docs.ext.intersphinx',
    'docs.ext.todo',
    'docs.ext.coverage',
    'docs.ext.imgmath',
    'docs.ext.mathjax',
    'docs.ext.ifconfig',
    'docs.ext.viewcode',
    'docs.ext.githubpages',
    'docs.ext.napoleon',
    'sphinx_rtd_theme',
    'docs.ext.autosummary',
    'sphinx_automodapi.automodapi'
]
numpydoc_show_class_members = False
graphviz_dot = "C:\\Users\\James\\Anaconda3\\Library\\bin\\graphviz\\dot.exe"


# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']


# -- Extension configuration -------------------------------------------------

# -- Options for intersphinx extension ---------------------------------------

# Example configuration for intersphinx: refer to the Python standard library.
intersphinx_mapping = {'https://docs.python.org/3/': None}

# -- Options for todo extension ----------------------------------------------

# If true, `todo` and `todoList` produce output, else they produce nothing.
todo_include_todos = True