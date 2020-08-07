import os

import numpy as np
from obspy import read

from lsforce import LSData

# Relative tolerance for test, see:
# https://numpy.org/doc/stable/reference/generated/numpy.testing.assert_allclose.html
RTOL = 3e-7

# Get location of test data
script_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(script_dir, 'data')

# Read in saved input ObsPy Stream
st_in = read(os.path.join(data_dir, 'st_in.pkl'), format='PICKLE')


def test_lsdata_st_proc():

    print('Testing LSData Stream processing...')

    # Create LSData object
    data = LSData(st_in, source_lat=60.0273, source_lon=-153.0683)

    # Test processed Stream
    for tr, tr_baseline in zip(
        data.st_proc, read(os.path.join(data_dir, 'st_proc.pkl'), format='PICKLE')
    ):
        print(f'Testing {tr.id}...')
        np.testing.assert_allclose(
            tr.data, tr_baseline.data, rtol=RTOL,
        )
