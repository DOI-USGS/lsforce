import datetime
import sys
from os import chdir, getcwd, pardir, path

root = path.abspath(path.join(path.dirname(path.abspath(__file__)), pardir))
sys.path.insert(0, root)
from versioneer import get_version

# Get version (this method avoids importing the package, see GitHub comment link below)
# https://github.com/python-versioneer/python-versioneer/issues/185#issue-343550345
owd = getcwd()
chdir(root)
__version__ = get_version()
chdir(owd)

SHOW_PRIVATE = False  # Set to True to build docs for private functions, methods, etc.

project = 'lsforce'

author = 'Kate E. Allstadt and Liam Toney'

copyright = f'{datetime.date.today().year}, {author}'

version = __version__
release = __version__

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.intersphinx',
    'recommonmark',
    'sphinx.ext.viewcode',
    'sphinxcontrib.apidoc',
    'sphinx.ext.mathjax',
    'sphinx.ext.todo',
]

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

html_theme = 'sphinx_rtd_theme'

html_title = f'{project} documentation'

napoleon_numpy_docstring = False

master_doc = 'index'

autodoc_mock_imports = ['cartopy', 'matplotlib', 'numpy', 'obspy', 'scipy', 'xarray']

apidoc_module_dir = path.join(root, 'lsforce')
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
    'xarray': ('http://xarray.pydata.org/en/stable/', None),
    'matplotlib': ('https://matplotlib.org/', None),
}
