import subprocess
from filecmp import cmp
from pathlib import Path

from lsforce import axisem2cps

# Get location of reference model file
ak135f_mod_ref = Path(__file__).resolve().parent / 'data' / 'ak135f.mod'


def test_axisem2cps_script():
    p = subprocess.run(['axisem2cps'], capture_output=True, text=True)
    assert not p.stderr


def test_axisem2cps_ak135f():
    axisem2cps('ak135f')
    try:
        assert cmp('ak135f.mod', ak135f_mod_ref, shallow=False)
    except AssertionError:
        Path('ak135f.mod').unlink()
        raise
    Path('ak135f.mod').unlink()
