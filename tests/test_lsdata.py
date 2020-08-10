import os
import pickle

import numpy as np
import pytest
from obspy import read

from lsforce import LSData

# Relative tolerance for test, see:
# https://numpy.org/doc/stable/reference/generated/numpy.testing.assert_allclose.html
RTOL = 3e-7

# Get location of test data
script_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(script_dir, 'data')

# Read in saved input ObsPy Stream
st_in = read(os.path.join(data_dir, 'input_stream.pkl'), format='PICKLE')

# Read in saved LSData object to test its plotting methods (note that this LSData object
# has its "st_orig" attribute removed to save space!)
with open(os.path.join(data_dir, 'lsdata.pkl'), 'rb') as f:
    lsdata = pickle.load(f)


def test_lsdata_st_proc():

    print('Testing LSData Stream processing...')

    # Create LSData object
    data = LSData(st_in, source_lat=60.0273, source_lon=-153.0683)

    # Test processed Stream
    for tr, tr_baseline in zip(
        data.st_proc,
        read(os.path.join(data_dir, 'processed_stream.pkl'), format='PICKLE'),
    ):
        print(f'Testing {tr.id}...')
        np.testing.assert_allclose(
            tr.data, tr_baseline.data, rtol=RTOL,
        )


@pytest.mark.mpl_image_compare(savefig_kwargs=dict(bbox_inches='tight'))
def test_lsdata_plot_data():
    return lsdata.plot_data(period_range=(15, 80))


@pytest.mark.mpl_image_compare(savefig_kwargs=dict(bbox_inches='tight'))
def test_lsdata_plot_stations():
    return lsdata.plot_stations(label_stations=True)
