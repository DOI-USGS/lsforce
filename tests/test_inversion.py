from lsforce import LSData, LSForce
from obspy import read
import numpy as np
import os
import shutil
import time

# Relative tolerance for test, see:
# https://numpy.org/doc/stable/reference/generated/numpy.testing.assert_allclose.html
RTOL = 3e-7

# Get start time
start_time = time.time()

# The directory where this script is located
script_dir = os.path.dirname(os.path.realpath(__file__))

# The directory where test data is located
data_dir = os.path.join(script_dir, 'data', 'test_inversion')

# Create a temporary run directory
tmp_dir = os.path.join(script_dir, 'tmp')
os.mkdir(tmp_dir)

# Read in Stream
st = read(os.path.join(data_dir, 'data.pkl'), format='PICKLE')

# Create LSData object
print('Creating LSData object...')
data = LSData(st, source_lat=60.0273, source_lon=-153.0683)
print('Done')

# Create LSForce object
force = LSForce(data=data, sampling_rate=1, nickname='test', main_folder=tmp_dir)

# Compute GFs
print('Computing Green\'s functions...')
force.compute_greens(
    model_file=os.path.join(data_dir, 'tak135sph.mod'), gf_duration=200, T0=-10
)
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
try:
    np.testing.assert_allclose(
        force.model,
        np.load(os.path.join(data_dir, 'model.npy')),
        verbose=False,
        rtol=RTOL,
    )
    print('PASS')
except AssertionError as error:
    print('FAIL')
    print(error.__str__() + '\n')

# Clean up
shutil.rmtree(tmp_dir)

# Print elapsed time (in seconds)
print(f'Elapsed time: {time.time() - start_time:.1f} seconds')
