import subprocess
from filecmp import cmp
from pathlib import Path

from lsforce import axisem2cps

# Get location of reference model file
iasp91_mod_ref = Path(__file__).resolve().parent / 'data' / 'iasp91.mod'


def test_axisem2cps_script():
    p = subprocess.run(['axisem2cps'], capture_output=True, text=True)
    assert not p.stderr


def test_axisem2cps_iasp91():
    axisem2cps('iasp91')
    try:
        assert cmp('iasp91.mod', iasp91_mod_ref, shallow=False)
    except AssertionError:
        Path('iasp91.mod').unlink()
        raise
    Path('iasp91.mod').unlink()
