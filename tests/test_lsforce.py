import copy
import os
import pickle
import tempfile

import numpy as np
import pytest
from obspy import read

from lsforce import LSForce, readrun
from lsforce.lsforce import GF_STARTTIME

# Set kwargs for all pytest-mpl tests
PYTEST_MPL_KWARGS = dict(style='default', savefig_kwargs=dict(bbox_inches='tight'))

# Toggle creation of test data
GENERATE_LSFORCE_PICKLE = False

# Relative tolerance for test, see:
# https://numpy.org/doc/stable/reference/generated/numpy.testing.assert_allclose.html
RTOL = 1e-6

# [m] Absolute tolerance for Green's function comparison
GF_ATOL = 1e-18

# Shared parameters for both inversion methods
SETUP_KWARGS = dict(period_range=(15, 80), syngine_model='iasp91_2s')
INVERT_KWARGS = dict(
    zero_time=119,
    impose_zero_start=True,
    add_to_zero=True,
    zero_start_taper_length=20,
    alpha=4.8e-17,
    zero_scaler=15,
    tikhonov_ratios=[0.4, 0.0, 0.6],
)

# Get location of test data
script_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(script_dir, 'data')

# Grab LSForce object from test data dir to use for LSForce plotting and forward model
# tests (and to grab its "data" attribute, which is an LSData object, to use as input
# for LSForce creation tests
lsforce = readrun(os.path.join(data_dir, 'lsforce.pkl'))


def test_lsforce_full():
    print('Testing full method...')

    # A temporary directory to run tests in
    temp_dir = tempfile.TemporaryDirectory()

    # Create LSForce object
    force = LSForce(
        data=lsforce.data,
        data_sampling_rate=1,
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

    # Generate a new LSForce object to use for testing data if requested
    if GENERATE_LSFORCE_PICKLE:
        with open(os.path.join(data_dir, 'lsforce_NEW.pkl'), 'wb') as f:
            pickle.dump(force, f)

    # Test resulting model
    print('Testing...')
    np.testing.assert_allclose(
        force.model, np.load(os.path.join(data_dir, 'model_full.npy')), rtol=RTOL
    )


def test_lsforce_triangle():
    print('Testing triangle method...')

    # A temporary directory to run tests in
    temp_dir = tempfile.TemporaryDirectory()

    # Create LSForce object
    force = LSForce(
        data=lsforce.data,
        data_sampling_rate=1,
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
        force.model, np.load(os.path.join(data_dir, 'model_triangle.npy')), rtol=RTOL
    )


def test_lsforce_gfs():
    print('Testing Syngine and CPS Green\'s functions...')

    # Create a single-station LSData object
    gf_data = copy.deepcopy(lsforce.data)
    for tr in gf_data.st_proc:
        if tr.stats.station != 'KALN':
            gf_data.st_proc.remove(tr)

    # Temporary directories to run tests in
    temp_dir_syn = tempfile.TemporaryDirectory()
    temp_dir_cps = tempfile.TemporaryDirectory()

    # Create LSForce objects
    force_syn = LSForce(
        data=gf_data, data_sampling_rate=1, main_folder=temp_dir_syn.name
    )
    force_cps = LSForce(
        data=gf_data, data_sampling_rate=1, main_folder=temp_dir_cps.name
    )

    # Different kwargs for CPS GFs
    setup_kwargs = copy.deepcopy(SETUP_KWARGS)
    del setup_kwargs['syngine_model']

    # Obtain GFs
    print('Grabbing GFs...')
    force_syn.setup(**SETUP_KWARGS)
    st_syn = force_syn.filtered_gf_st
    force_cps.setup(cps_model=os.path.join(data_dir, 'iasp91.mod'), **setup_kwargs)
    st_cps = force_cps.filtered_gf_st
    print('Done')

    # Clean up temporary dirs
    temp_dir_syn.cleanup()
    temp_dir_cps.cleanup()

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


def test_forward():
    print('Testing forward problem...')

    st_syn = lsforce.forward(lsforce.Z, lsforce.N, lsforce.E)

    # Test that st_syn data are IDENTICAL to dtnew
    np.testing.assert_allclose(
        np.array([tr.data for tr in st_syn]),
        lsforce.dtnew,
    )


pytest_mpl_kwargs = copy.deepcopy(PYTEST_MPL_KWARGS)
pytest_mpl_kwargs['savefig_kwargs']['bbox_inches'] = None  # Plot already sized nicely
pytest_mpl_kwargs['savefig_kwargs']['dpi'] = 200  # Since this plot is small


@pytest.mark.mpl_image_compare(**pytest_mpl_kwargs)
def test_lsforce_plot_fits():
    return lsforce.plot_fits(xlim=(-50, 250))


@pytest.mark.mpl_image_compare(**PYTEST_MPL_KWARGS)
def test_lsforce_plot_fits_legacy():
    return lsforce.plot_fits(xlim=(-50, 250), legacy=True)


@pytest.mark.mpl_image_compare(**PYTEST_MPL_KWARGS)
def test_lsforce_plot_forces():
    return lsforce.plot_forces(xlim=(-50, 250))


@pytest.mark.mpl_image_compare(**PYTEST_MPL_KWARGS)
def test_lsforce_plot_angle_magnitude():
    return lsforce.plot_angle_magnitude(xlim=(-50, 250))
