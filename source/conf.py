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


# -- Project information -----------------------------------------------------
master_doc = 'index'

project = "Single cell analysis"
copyright = "2022, Cameron Smith"
author = "Cameron Smith"


# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named "sphinx.ext.*") or your custom
# ones.
extensions = [
    "nbsphinx",
    "sphinx.ext.mathjax",
    "sphinxcontrib.bibtex",        # for bibliographic references
    "sphinx_copybutton",           # for adding “copy to clipboard” buttons to all text/code boxes | commented due to multiple scrollbar issue https://github.com/cameronraysmith/sc/issues/1
    "sphinxcontrib.rsvgconverter", # for SVG->PDF conversion in LaTeX output
    "sphinx_design",
    "sphinx_inline_tabs",
]

# Add any paths that contain templates here, relative to this directory.
# templates_path = ["_templates"]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['*.py','*.md', '*/*.py','*/*.md']
nbsphinx_execute = "never"

pygments_style = "sphinx"
pygments_dark_style = "monokai"

html_theme_options = {
    "footer_icons": [
        {
            "name": "GitHub",
            "url": "https://github.com/cameronraysmith/sc",
            "html": """
                <svg stroke="currentColor" fill="currentColor" stroke-width="0" viewBox="0 0 16 16">
                    <path fill-rule="evenodd" d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0 0 16 8c0-4.42-3.58-8-8-8z"></path>
                </svg>
            """,
            "class": "",
        },
    ],
    "source_repository": "https://github.com/cameronraysmith/sc/",
    "source_branch": "master",
    "source_directory": "source/",
}

# <script>
#         if (document.location.host) {
#           $(document.currentScript).replaceWith(
#             '<a class="reference external" ' +
#             'href="https://nbviewer.jupyter.org/url' +
#             (window.location.protocol == 'https:' ? 's/' : '/') +
#             window.location.host +
#             window.location.pathname.slice(0, -4) +
#             'ipynb"><em>nbviewer</em></a>.'
#           );
#         }
# </script>

# This is processed by Jinja2 and inserted before each notebook
nbsphinx_prolog = r"""
{% set docname = 'source/' + env.doc2path(env.docname, base=None) %}

.. raw:: html

    <div class="admonition note">
      This page was generated from
      <a class="reference external" href="https://github.com/cameronraysmith/sc/blob/master/{{ docname|e }}">{{ docname|e }}</a>.
      <span style="white-space: nowrap;">
      <a href="https://mybinder.org/v2/gh/cameronraysmith/sc/master?filepath={{ docname|e }}"><img alt="Binder badge" src="https://mybinder.org/badge_logo.svg" style="vertical-align:text-bottom"></a> or 
      <a href="https://colab.research.google.com/github/cameronraysmith/sc/blob/master/{{ docname|e }}"><img alt="Colab badge" src="https://colab.research.google.com/assets/colab-badge.svg" style="vertical-align:text-bottom"></a>
      </span>
    </div>

.. raw:: latex

    \nbsphinxstartnotebook{\scriptsize\noindent\strut
    \textcolor{gray}{The following section was generated from
    \sphinxcode{\sphinxupquote{\strut {{ docname | escape_latex }}}} \dotfill}}
"""

# This is processed by Jinja2 and inserted after each notebook
nbsphinx_epilog = r"""
{% set docname = 'source/' + env.doc2path(env.docname, base=None) %}
.. raw:: latex

    \nbsphinxstopnotebook{\scriptsize\noindent\strut
    \textcolor{gray}{\dotfill\ \sphinxcode{\sphinxupquote{\strut
    {{ docname | escape_latex }}}} ends here.}}
"""



mathjax3_config = {
    'tex': {'tags': 'ams', 'useLabelIds': True},
}

bibtex_bibfiles = ['references.bib']

# Support for notebook formats other than .ipynb
nbsphinx_custom_formats = {
    '.pct.py': ['jupytext.reads', {'fmt': 'py:percent'}],
    '.md': ['jupytext.reads', {'fmt': 'md'}],
}

# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
# html_theme = 'alabaster'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
# html_static_path = ['_static']

# -- Get version information and date from Git ----------------------------

try:
    from subprocess import check_output
    release = check_output(['git', 'describe', '--tags', '--always'])
    release = release.decode().strip()
    today = check_output(['git', 'show', '-s', '--format=%ad', '--date=short'])
    today = today.decode().strip()
except Exception:
    release = ''
    today = ''

# -- Options for HTML output ----------------------------------------------

# html_favicon = 'favicon.svg'
# html_title = project + ' version ' + release
html_title = project

# -- Options for LaTeX output ---------------------------------------------

# See https://www.sphinx-doc.org/en/master/latex.html
latex_elements = {
    'papersize': 'a4paper',
    'printindex': '',
    'sphinxsetup': r"""
        %verbatimwithframe=false,
        %verbatimwrapslines=false,
        %verbatimhintsturnover=false,
        VerbatimColor={HTML}{F5F5F5},
        VerbatimBorderColor={HTML}{E0E0E0},
        noteBorderColor={HTML}{E0E0E0},
        noteborder=1.5pt,
        warningBorderColor={HTML}{E0E0E0},
        warningborder=1.5pt,
        warningBgColor={HTML}{FBFBFB},
    """,
    'preamble': r"""
\usepackage[sc,osf]{mathpazo}
\linespread{1.05}  % see http://www.tug.dk/FontCatalogue/urwpalladio/
\renewcommand{\sfdefault}{pplj}  % Palatino instead of sans serif
\IfFileExists{zlmtt.sty}{
    \usepackage[light,scaled=1.05]{zlmtt}  % light typewriter font from lmodern
}{
    \renewcommand{\ttdefault}{lmtt}  % typewriter font from lmodern
}
\usepackage{booktabs}  % for Pandas dataframes
""",
}

latex_documents = [
    (master_doc, 'nbsphinx.tex', project, author, 'howto'),
]

latex_show_urls = 'footnote'
latex_show_pagerefs = True

# -- Options for EPUB output ----------------------------------------------

# These are just defined to avoid Sphinx warnings related to EPUB:
version = release
suppress_warnings = ['epub.unknown_project_files']

# -- Set default HTML theme (if none was given above) ---------------------
# html_theme = 'sphinx_book_theme'

if 'html_theme' not in globals():
    try:
        # import insipid_sphinx_theme
        import furo
    except ImportError:
        pass
    else:
        # html_theme = 'insipid'
        html_theme = 'furo'
        # html_copy_source = False
        html_permalinks_icon = '\N{SECTION SIGN}'