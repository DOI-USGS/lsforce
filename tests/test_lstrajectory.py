import os

import pytest

from lsforce import LSTrajectory, readrun

# Set kwargs for all pytest-mpl tests
PYTEST_MPL_KWARGS = dict(style='default', savefig_kwargs=dict(bbox_inches='tight'))

# Get location of test data
script_dir = os.path.dirname(os.path.realpath(__file__))
data_dir = os.path.join(script_dir, 'data')

# Grab LSForce object from test data dir to use for LSTrajectory creation
lsforce = readrun(os.path.join(data_dir, 'lsforce.pkl'))


@pytest.mark.mpl_image_compare(**PYTEST_MPL_KWARGS)
def test_lstrajectory_calculation_and_plots():

    print('Testing LSTrajectory calculation and plots...')

    # Create LSTrajectory object (which also calculates trajectory)
    trajectory = LSTrajectory(
        lsforce, target_length=6, duration=200, detrend_velocity=200
    )

    # Return map-view figure
    return trajectory.plot_trajectory()
