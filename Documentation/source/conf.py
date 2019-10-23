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
sys.path.insert(0, os.path.abspath('../..'))
sys.setrecursionlimit(1500)


# -- Project information -----------------------------------------------------

project = 'polymermassspec'
copyright = '2019, David Doran, Graham Keenan, Emma Clarke, Cole Mathis and Leroy Cronin'
author = 'David Doran, Graham Keenan, Emma Clarke, Cole Mathis and Leroy Cronin'

# The full version, including alpha/beta/rc tags
release = '0.0.1'


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = ['sphinx.ext.autodoc', 'sphinx.ext.coverage', 'sphinx.ext.napoleon'
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'classic'

html_theme_options = {

'stickysidebar': 'true',

# Makes sidebar collapsable.
'collapsiblesidebar':'false',

# Background color for the footer line.
'footerbgcolor':'White',

# Text color for the footer line.
'footertextcolor':'Black',

# Background color for the sidebar.
'sidebarbgcolor':'LightGray',

# Background color for the sidebar collapse button.
'sidebarbtncolor':'DarkGrey', 

# Text color for the sidebar 
'sidebartextcolor':'Black', 

# Link color for the sidebar
'sidebarlinkcolor':'Black',

# Background color for the relation bar.
'relbarbgcolor':'LightBlue',

# Text color for the relation bar.
'relbartextcolor':'White',

# Body background color.
'bgcolor':'WhiteSmoke',

# Body text color.
'textcolor':'Black',

# Body link color
'linkcolor':'SteelBlue',

# Body color for visited links.
'visitedlinkcolor':'SteelBlue',

# Background color for headings.
'headbgcolor':'WhiteSmoke',

# Text color for headings
'headtextcolor':'Black',

# Link color for headings.
'headlinkcolor':'Black',

# Background color for code blocks.
'codebgcolor':'WhiteSmoke',

# Default text color for code blocks, if not set differently by the highlighting style.
'codetextcolor':'Black',

# Font for normal text.
'bodyfont':'Arial',

# Font for headings.
'headfont':'Arial'
 }

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']