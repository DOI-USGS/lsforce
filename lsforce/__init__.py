from importlib.metadata import version

from .lsdata import LSData, make_lsdata_syn
from .lsforce import LSForce, readrun
from .lstrajectory import LSTrajectory

__version__ = version('lsforce')
del version
