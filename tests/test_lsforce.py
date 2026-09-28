import copy
import os
import pickle
import tempfile

import numpy as np
import pandas as pd
import pytest

from lsforce import LSForce, readrun
from lsforce.lsforce import GF_STARTTIME

# Set kwargs for all pytest-mpl tests
PYTEST_MPL_KWARGS = dict(style='default', savefig_kwargs=dict(bbox_inches='tight'))

# Toggle creation of test data
GENERATE_LSFORCE_PICKLE = False

# Relative tolerance for test, see:
# https://numpy.org/doc/stable/reference/generated/numpy.testing.assert_allclose.html
RTOL = 1e-5

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


def test_lsforce_jackknifing():
    print('Testing jackknifing...')

    # Add fake additional channels to LSData object to test jackknifing
    data = copy.deepcopy(lsforce.data)
    for station in set(tr.stats.station for tr in data.st_proc):
        st_sta = data.st_proc.select(station=station)
        for component in 'ZRT':
            if not st_sta.select(component=component):
                _tr = st_sta[0].copy()  # Template trace
                _tr.stats.component = component
                data.st_proc += _tr
    data.st_proc.sort(keys=['distance', 'station', 'channel'])

    # A temporary directory to run tests in
    temp_dir = tempfile.TemporaryDirectory()

    # Create LSForce object
    force = LSForce(
        data=data,
        data_sampling_rate=1,
        main_folder=temp_dir.name,
        method='full',
    )

    # Needed for some tests below
    unique_station_list = pd.unique(  # Preserve order of stations
        np.array([tr.stats.station for tr in force.data.st_proc])
    ).tolist()

    # Set up for inversion
    print('Setting up...')
    force.setup(**SETUP_KWARGS)
    print('Done')

    # (A) Channel-based leave-one-out analysis
    force.invert(jackknife='channel', frac_delete=1, **INVERT_KWARGS)
    assert force.jackknife.num_iter == force.data.st_proc.count()
    all_tr_ids = [tr.id for tr in force.data.st_proc]
    missing_tr_ids = []
    for tr_id_list in force.jackknife.tr_id_lists_all:
        for tr_id in all_tr_ids:
            if tr_id not in tr_id_list:
                missing_tr_ids.append(tr_id)
    assert missing_tr_ids == all_tr_ids  # Should match, including order!

    # (B) Station-based leave-one-out analysis
    force.invert(jackknife='station', frac_delete=1, **INVERT_KWARGS)
    assert force.jackknife.num_iter == len(unique_station_list)
    missing_stations = []
    for tr_id_list in force.jackknife.tr_id_lists_all:
        stations = []
        for tr_id in tr_id_list:
            stations.append(tr_id.split('.')[1])  # Get station from trace ID
        for station in unique_station_list:
            if station not in set(stations):
                missing_stations.append(station)
    assert missing_stations == unique_station_list  # Should match, including order!

    # (C) Channel-based jackknife iterations
    force.invert(jackknife='channel', num_iter=10, frac_delete=0.5, **INVERT_KWARGS)
    expected_n_chans = force.data.st_proc.count() - int(
        round(force.jackknife.frac_delete * force.data.st_proc.count())
    )
    n_chans = set([len(tr_id_list) for tr_id_list in force.jackknife.tr_id_lists_all])
    assert (len(n_chans) == 1) and (list(n_chans)[0] == expected_n_chans)

    # (D) Station-based jackknife iterations
    force.invert(jackknife='station', num_iter=10, frac_delete=0.5, **INVERT_KWARGS)
    expected_n_stations = len(unique_station_list) - int(
        round(force.jackknife.frac_delete * len(unique_station_list))
    )
    n_stations = []
    for tr_id_list in force.jackknife.tr_id_lists_all:
        stations = []
        for tr_id in tr_id_list:
            stations.append(tr_id.split('.')[1])  # Get station from trace ID
        n_stations.append(len(set(stations)))
    n_stations = set(n_stations)
    assert (len(n_stations) == 1) and (list(n_stations)[0] == expected_n_stations)

    # Clean up temporary dir
    temp_dir.cleanup()


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
