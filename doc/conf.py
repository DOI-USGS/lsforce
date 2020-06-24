import os
import sys

SHOW_PRIVATE = False  # Set to True to build docs for private functions, methods, etc.

sys.path.insert(0, os.path.abspath('..'))

project = 'lsforce'

html_show_copyright = False

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.intersphinx',
    'recommonmark',
    'sphinx.ext.viewcode',
    'sphinxcontrib.apidoc',
    'sphinx.ext.mathjax',
]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

html_theme = 'sphinx_rtd_theme'

napoleon_numpy_docstring = False

master_doc = 'index'

apidoc_module_dir = '../lsforce'
apidoc_output_dir = 'api'
apidoc_separate_modules = True
apidoc_toc_file = False

autodoc_default_options = {
    'special-members': '__init__',
    'undoc-members': True,
}
if SHOW_PRIVATE:
    autodoc_default_options['private-members'] = True

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'obspy': ('https://docs.obspy.org/', None),
    'xarray': ('http://xarray.pydata.org/en/stable/', None),
    'matplotlib': ('https://matplotlib.org/', None),
}
