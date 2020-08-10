import copy
import os
import pickle
import tempfile

import numpy as np
from obspy import read

from lsforce import LSForce
from lsforce.lsforce import GF_STARTTIME

# Relative tolerance for test, see:
# https://numpy.org/doc/stable/reference/generated/numpy.testing.assert_allclose.html
RTOL = 1e-6

# [m] Absolute tolerance for Green's function comparison
GF_ATOL = 1e-18

# Shared parameters for both inversion methods
SETUP_KWARGS = dict(period_range=(15, 80), zerophase=True, syngine_model='iasp91_2s')
INVERT_KWARGS = dict(
    zero_time=119,
    impose_zero=True,
    add_to_zero=True,
    alpha=4.8e-17,
    zero_scaler=15,
    tikhonov_ratios=[0.4, 0.0, 0.6],
)

# Get location of test data
script_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(script_dir, 'data')

# Read in saved LSData object, since we're just testing LSForce here (note that this
# LSData object has it's "st_orig" attribute removed to save space!)
with open(os.path.join(data_dir, 'lsdata.pkl'), 'rb') as f:
    data = pickle.load(f)


def test_lsforce_full():

    print('Testing full method...')

    # A temporary directory to run tests in
    temp_dir = tempfile.TemporaryDirectory()

    # Create LSForce object
    force = LSForce(
        data=data,
        data_sampling_rate=1,
        nickname='full',
        main_folder=temp_dir.name,
        method='full',
    )

    # Set up for inversion
    print('Setting up...')
    force.setup(**SETUP_KWARGS)
    print('Done')

    # Invert
    print('Inverting...')
    force.invert(**INVERT_KWARGS)
    print('Done')

    # Clean up temporary dir
    temp_dir.cleanup()

    # Test resulting model
    print('Testing...')
    np.testing.assert_allclose(
        force.model, np.load(os.path.join(data_dir, 'model_full.npy')), rtol=RTOL,
    )


def test_lsforce_triangle():

    print('Testing triangle method...')

    # A temporary directory to run tests in
    temp_dir = tempfile.TemporaryDirectory()

    # Create LSForce object
    force = LSForce(
        data=data,
        data_sampling_rate=1,
        nickname='triangle',
        main_folder=temp_dir.name,
        method='triangle',
    )

    # Set up for inversion
    print('Setting up...')
    force.setup(triangle_half_width=10, **SETUP_KWARGS)
    print('Done')

    # Invert
    print('Inverting...')
    force.invert(**INVERT_KWARGS)
    print('Done')

    # Clean up temporary dir
    temp_dir.cleanup()

    # Test resulting model
    print('Testing...')
    np.testing.assert_allclose(
        force.model, np.load(os.path.join(data_dir, 'model_triangle.npy')), rtol=RTOL,
    )


def test_lsforce_gfs():

    print('Testing Syngine and CPS Green\'s functions...')

    # Create a single-station LSData object
    gf_data = copy.copy(data)
    for tr in gf_data.st_proc:
        if tr.stats.station != 'KALN':
            gf_data.st_proc.remove(tr)

    # A temporary directory to run tests in
    temp_dir = tempfile.TemporaryDirectory()

    # Create LSForce object
    force = LSForce(
        data=gf_data, data_sampling_rate=1, nickname='gfs', main_folder=temp_dir.name,
    )

    # Obtain GFs
    print('Grabbing GFs...')
    force.setup(**SETUP_KWARGS)
    st_syn = force.filtered_gf_st
    st_cps = read(
        os.path.join(data_dir, 'filtered_gf_stream.pkl')
    )  # Since no CPS install
    print('Done')

    # Clean up temporary dir
    temp_dir.cleanup()

    # Make sure they have equal number of samples, no padding here since that will throw
    # off the test
    syn_npts = st_syn[0].stats.npts
    cps_npts = st_cps[0].stats.npts
    if syn_npts > cps_npts:
        st_syn.trim(GF_STARTTIME, st_cps[0].stats.endtime)
    elif cps_npts > syn_npts:
        st_cps.trim(GF_STARTTIME, st_syn[0].stats.endtime)

    # Test GFs to make sure they're about the same
    for tr_syn, tr_cps in zip(st_syn, st_cps):
        print(f'Testing {tr_syn.id}...')
        np.testing.assert_allclose(tr_syn.data, tr_cps.data, atol=GF_ATOL)
