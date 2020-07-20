import os
import pickle
import tempfile

import numpy as np

from lsforce import LSForce

# Relative tolerance for test, see:
# https://numpy.org/doc/stable/reference/generated/numpy.testing.assert_allclose.html
RTOL = 3e-7

# Shared parameters for both inversion methods
SETUP_KWARGS = dict(
    period_range=(15, 80),
    zerophase=True,
    syngine_model='iasp91_2s',
    gf_duration=200,
    T0=-10,
)
INVERT_KWARGS = dict(
    zero_time=119,
    impose_zero=True,
    add_to_zero=True,
    alphaset=4.8e-17,
    zero_scaler=2,
    tikhonov_ratios=[0.4, 0.0, 0.6],
)

# Get location of test data
script_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(script_dir, 'data', 'test_lsforce')

# Read in saved LSData object, since we're just testing LSForce here (note that this
# LSData object has it's "st_orig" attribute removed to save space!)
with open(os.path.join(data_dir, 'lsdata.pkl'), 'rb') as f:
    data = pickle.load(f)


def test_lsforce_tik():

    print('Testing tik method...')

    # A temporary directory to run tests in
    temp_dir = tempfile.TemporaryDirectory()

    # Create LSForce object
    force = LSForce(
        data=data,
        data_sampling_rate=1,
        nickname='tik',
        main_folder=temp_dir.name,
        method='tik',
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
        force.model, np.load(os.path.join(data_dir, 'model_tik.npy')), rtol=RTOL,
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
