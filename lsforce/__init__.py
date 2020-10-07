from ._version import get_versions
from .lsdata import LSData
from .lsforce import LSForce, readrun
from .lstrajectory import LSTrajectory

__version__ = get_versions()['version']
del get_versions
