# -*- coding: utf-8 -*-
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

sys.path.insert(0, os.path.abspath("../.."))


# -- Project information -----------------------------------------------------

project = "mofchecker"
copyright = "2022, Kevin Maik Jablonka"
author = "Kevin Maik Jablonka"

# The full version, including alpha/beta/rc tags
release = "0.9.6-dev"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.coverage",
    "sphinx.ext.napoleon",
    "sphinx_autodoc_typehints",
    "sphinx.ext.mathjax",
    "sphinx_copybutton",
    "sphinx.ext.viewcode",
    "sphinx.ext.autosectionlabel",
]

copybutton_selector = "div:not(.no-copy)>div.highlight pre"
copybutton_prompt_text = r">>> |\.\.\. |\$ |In \[\d*\]: | {2,5}\.\.\.: | {5,8}: "
copybutton_prompt_is_regexp = True

autosummary_generate = True

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


numfig = True


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "furo"
if os.path.exists("logo.png"):
    html_logo = "logo.png"

# Register the theme as an extension to generate a sitemap.xml


html_theme_options = {
    "path_to_docs": "docs",
    "repository_url": "https://github.com/kjappelbaum/mofchecker",
    "repo_url": "https://github.com/kjappelbaum/mofchecker",
    "repo_name": "mofchecker",
    "repo_type": "github",
    # "source_repository": "https://github.com/pradyunsg/furo/",
    "source_branch": "main",
    "source_directory": "docs/",
    "edit_uri": "https://github.com/kjappelbaum/mofchecker/docs/source",
}


html_show_sphinx = False


# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "sphinx"

# If true, SmartyPants will be used to convert quotes and dashes to
# typographically correct entities.
html_use_smartypants = True
