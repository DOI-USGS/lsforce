import os

import numpy as np
from obspy import read

from lsforce import LSData, LSForce

# Relative tolerance for test, see:
# https://numpy.org/doc/stable/reference/generated/numpy.testing.assert_allclose.html
RTOL = 3e-7


def test_inversion():

    # The directory where this script is located
    script_dir = os.path.dirname(os.path.realpath(__file__))

    # The directory where test data is located
    data_dir = os.path.join(script_dir, 'data', 'test_inversion')

    # Read in Stream
    st = read(os.path.join(data_dir, 'data.pkl'), format='PICKLE')

    # Create LSData object
    print('Creating LSData object...')
    data = LSData(st, source_lat=60.0273, source_lon=-153.0683)
    print('Done')

    # Create LSForce object
    force = LSForce(
        data=data, data_sampling_rate=1, nickname='test', main_folder=data_dir
    )

    # Load GFs
    print('Loading Green\'s functions...')
    force.load_greens(model_file=os.path.join(data_dir, 'tak135sph.mod'))
    print('Done')

    # Set up for inversion
    force.setup(period_range=(15, 80), zerophase=True)

    # Invert
    print('Inverting...')
    force.invert(
        zero_time=119,
        impose_zero=True,
        add_to_zero=True,
        alphaset=4.8e-20,
        zero_scaler=2,
        tikhonov_ratios=[0.4, 0.0, 0.6],
    )
    print('Done')

    # Test resulting model
    print('Testing...')
    np.testing.assert_allclose(
        force.model, np.load(os.path.join(data_dir, 'model.npy')), rtol=RTOL,
    )
