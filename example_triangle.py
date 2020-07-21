import os

from obspy import UTCDateTime, read
from obspy.clients.fdsn import Client

from lsforce import LSData, LSForce

# Arbitrary run directory containing model file
LSFORCE_RUN_DIR = os.getcwd()

RUN_NAME = 'hunza_attabad'  # Nickname for this run

PERIOD_RANGE = (50, 150)  # [s] Bandpass filter corners

LS_LAT, LS_LON = (36.31, 74.82)  # Where the point force will be applied
ORIGIN_TIME = UTCDateTime(2010, 1, 4, 8, 36)  # From Ekström supp table

STARTTIME = ORIGIN_TIME - 100
ENDTIME = ORIGIN_TIME + 600

# Set up folder structure
main_folder = os.path.join(LSFORCE_RUN_DIR, RUN_NAME)
if not os.path.exists(main_folder):
    os.mkdir(main_folder)

# If GF's don't exist yet, or if the Stream has been modified, this must be True
CALCULATE_GF = True

#%% GATHER INVERSION WAVEFORMS

data_filename = os.path.join(main_folder, f'{RUN_NAME}_data.pkl')

# Download data if it doesn't exist as a file
if not os.path.exists(data_filename):

    client = Client('IRIS')
    waveform_kwargs = dict(
        channel='BHE,BHN,BHZ',
        location='00,',
        starttime=STARTTIME,
        endtime=ENDTIME,
        attach_response=True,
    )

    NETWORKS = (
        'II',
        'IU',
        'KN',
    )
    STATIONS = (
        'KBL',
        'KZA',
        'AAK',
        'EKS2',
    )
    st = client.get_waveforms(
        network=','.join(NETWORKS), station=','.join(STATIONS), **waveform_kwargs,
    )

    # Remove extra channels for station AAK
    for tr in st.select(station='AAK', location=''):
        st.remove(tr)

    # Remove horizontals
    for station in ['KZA', 'EKS2']:
        for tr in st.select(station=station, component='[EN]'):
            st.remove(tr)

    # Grab coordinates
    inv = client.get_stations(
        network=','.join(NETWORKS),
        starttime=STARTTIME,
        endtime=ENDTIME,
        level='channel',
    )

    # Assign coordinates to Traces
    for tr in st:
        coords = inv.get_coordinates(tr.id, datetime=STARTTIME)
        tr.stats.latitude = coords['latitude']
        tr.stats.longitude = coords['longitude']

    st.write(data_filename, format='PICKLE')

# Use file if it exists, for speed
else:
    st = read(data_filename, format='PICKLE')

# Create LSData object
data = LSData(st, source_lat=LS_LAT, source_lon=LS_LON)

# Detrend and taper intensively!
data.st_proc.detrend('polynomial', order=12)
data.st_proc.taper(0.3)

# Now that it's rotated, remove some radial components
for station in ['ABKT', 'KIV', 'ANTO']:
    for tr in data.st_proc.select(station=station, component='R'):
        data.st_proc.remove(tr)

# THESE STATIONS WERE INCLUDED IN EKSTRÖM BUT LOOK BAD TO ME
for tr in data.st_proc.select(station='ANTO', component='Z'):
    data.st_proc.remove(tr)
for tr in data.st_proc.select(station='ABKT', component='T'):
    data.st_proc.remove(tr)
for tr in data.st_proc.select(station='KIV', component='T'):
    data.st_proc.remove(tr)

# Create plots
data.plot_stations(label_stations=True, gshhs_scale='intermediate')
data.plot_data(equal_scale=True, period_range=PERIOD_RANGE)

#%% SETUP

force = LSForce(
    data=data,
    data_sampling_rate=1,
    nickname=RUN_NAME,
    main_folder=main_folder,
    method='triangle',
)

force.setup(
    period_range=PERIOD_RANGE,
    zerophase=True,
    syngine_model='iasp91_2s',
    triangle_half_width=10,
)

#%% INVERT

force.invert(
    zero_time=150,
    max_duration=65,
    impose_zero=True,
    add_to_zero=True,
    alphaset=3.6e-18,
    zero_scaler=1,
    tikhonov_ratios=(0, 0, 1),
)

#%% PLOT INVERSION

# Plot waveform fits
force.plot_fits()

# Plot results
force.plot_forces(xlim=[-50, 150])
