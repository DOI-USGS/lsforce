from datetime import datetime
from importlib.metadata import version

SHOW_PRIVATE = False  # Set to True to build docs for private functions, methods, etc.

project = 'lsforce'

author = 'Kate E. Allstadt and Liam Toney'

copyright = f'{datetime.now().year}, {author}'

__version__ = version('lsforce')
version = __version__
release = __version__

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.intersphinx',
    'recommonmark',
    'sphinx.ext.viewcode',
    'sphinxcontrib.apidoc',
    'sphinx.ext.todo',
    'sphinx_markdown_builder',
]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

html_theme = 'sphinx_rtd_theme'

html_title = f'{project} documentation'

html_theme_options = dict(display_version=False)

napoleon_numpy_docstring = False

root_doc = 'index'

apidoc_module_dir = '../lsforce'
apidoc_output_dir = 'api'
apidoc_separate_modules = True
apidoc_toc_file = False
apidoc_module_first = True

autoclass_content = 'both'

autodoc_default_options = {
    'undoc-members': True,
}
if SHOW_PRIVATE:
    autodoc_default_options['private-members'] = True

todo_include_todos = True

intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'obspy': ('https://docs.obspy.org/', None),
    'xarray': ('https://docs.xarray.dev/en/stable/', None),
    'matplotlib': ('https://matplotlib.org/', None),
}
