# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'wham_derivatives'
copyright = '2022, Zeke Piskulich'
author = 'Zeke Piskulich'
release = 'v0.9dev'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

import os, sys
sys.path.insert(0, os.path.abspath('../../wham_derivatives'))
extensions = ['sphinx.ext.napoleon','sphinx.ext.autodoc','sphinx.ext.todo','sphinx.ext.githubpages']

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']
