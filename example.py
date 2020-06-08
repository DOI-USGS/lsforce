from lsforce import LSforce
from waveform_collection import gather_waveforms
from obspy import read, UTCDateTime
from obspy.geodetics import gps2dist_azimuth
import os
import numpy as np

RUN_NAME = 'iliamna_2016'  # Nickname for this run
EVENT_ID = 1  # Provide this number so database-based code doesn't break

PERIOD_RANGE = (15, 80)  # [s] Bandpass filter corners

LS_LAT, LS_LON = (60.0273, -153.0683)  # Where the point force will be applied
ORIGIN_TIME = UTCDateTime(2016, 5, 22, 7, 57, 34)

STARTTIME = ORIGIN_TIME - 100
ENDTIME = ORIGIN_TIME + 300

MODEL_FILE = os.path.join(os.getcwd(), 'meow', 'tak135sph.mod')
MAIN_FOLDER = os.path.join(os.getcwd(), 'meow', RUN_NAME)
AUX_FOLDER = os.path.join(MAIN_FOLDER, f'EV{EVENT_ID}')

# Make folders if they don't already exist
if not os.path.exists(MAIN_FOLDER):
    os.mkdir(MAIN_FOLDER)
if not os.path.exists(AUX_FOLDER):
    os.mkdir(AUX_FOLDER)

# If GF's don't exist yet, or if the Stream has been modified, this must be True
CALCULATE_GF = True

#%% GATHER INVERSION WAVEFORMS

TONEY_ET_AL_NUM_CHANS = 32  # Number of channels used in paper's 2016 Iliamna inversion

data_filename = os.path.join(AUX_FOLDER, f'{RUN_NAME}_data.pkl')

# Download data if it doesn't exist as a file
if not os.path.exists(data_filename):

    # Gather vertical components (most of the waveforms!)
    NETWORKS = (
        'AK', 'AT', 'AV', 'TA', 'ZE',
    )
    STATIONS = (
        'KALN', 'HOM',  'HLC5', 'CLAM', 'SALA', 'LTUY', 'N19K', 'CNP',  'NSKI', 'LTUX',
        'BRLK', 'Q19K', 'LTUW', 'BRSE', 'WFLS', 'P18K', 'BING', 'CONG', 'WFLW', 'MPEN',
        'BULG', 'SLK',  'N18K', 'JOES', 'SVW2', 'JUDD', 'O22K', 'FIRE',
    )
    st = gather_waveforms(
        source='IRIS',
        network=','.join(NETWORKS),
        station=','.join(STATIONS),
        location='*',
        channel='BHZ,HHZ',
        starttime=STARTTIME,
        endtime=ENDTIME,
    )

    # Gather horizontals (only a few)
    st += gather_waveforms(
        source='IRIS',
        network='TA',
        station='N19K',
        location='*',
        channel='BHE,BHN',
        starttime=STARTTIME,
        endtime=ENDTIME,
    )
    st += gather_waveforms(
        source='IRIS',
        network='ZE',
        station='WFLW',
        location='*',
        channel='HHE,HHN',
        starttime=STARTTIME,
        endtime=ENDTIME,
    )

    # Assign additional info to Traces
    for tr in st:
        dist, az, baz = gps2dist_azimuth(
            LS_LAT, LS_LON, tr.stats.latitude, tr.stats.longitude
        )
        tr.stats.rdist = dist / 1000  # [km]
        tr.stats.azimuth = az  # [deg]
        tr.stats.back_azimuth = baz  # [deg]

    # Rotate (on a station-by-station basis!)
    for station in np.unique([tr.stats.station for tr in st]):
        st.select(station=station).rotate('NE->RT')

    st.detrend('polynomial', order=2)
    st.remove_response(output='DISP', water_level=60, zero_mean=False)

    st.write(data_filename, format='PICKLE')

# Use file if it exists, for speed
else:
    st = read(data_filename, format='PICKLE')

# Verify that the correct number of channels has been retrieved from IRIS
assert st.count() == TONEY_ET_AL_NUM_CHANS, 'Not the correct number of channels.'

st.sort(keys=['rdist', 'channel'])

#%% GATHER REFERENCE WAVEFORMS

RAYLEIGH_VELO = 0.9  # [km/s] Surface-wave group velocity @ 1 Hz
INFRA_VELO = 0.337   # [km/s] Reasonable given air temp of 50 degrees F

# Gather seismic
st_hf = gather_waveforms(
    source='IRIS',
    network='AV',
    station='ILSW',
    location='--',
    channel='BHZ',
    starttime=STARTTIME,
    endtime=ENDTIME,
)
# Gather infrasound
st_infra = gather_waveforms(
    source='IRIS',
    network='TA',
    station='O20K',
    location='*',
    channel='BDF',
    starttime=STARTTIME,
    endtime=ENDTIME,
)

# Combined processing
(st_hf + st_infra).remove_response()
(st_hf + st_infra).detrend('demean')
(st_hf + st_infra).taper(max_percentage=0.05)

# Separate filtering
st_hf.filter('bandpass', freqmin=0.5, freqmax=5)
st_infra.filter('bandpass', freqmin=0.5, freqmax=10)

# Add "rdist" to tr.stats
for tr in st_hf + st_infra:
    dist = gps2dist_azimuth(LS_LAT, LS_LON, tr.stats.latitude, tr.stats.longitude)[0]
    tr.stats.rdist = dist / 1000  # [km]

# Approximate correction for travel time
hf_shift = st_hf[0].stats.rdist / RAYLEIGH_VELO
infra_shift = st_infra[0].stats.rdist / INFRA_VELO

#%% SETUP

force = LSforce(
    st=st,
    samplerate=1,
    nickname=RUN_NAME,
    event_id=EVENT_ID,
    mainfolder=MAIN_FOLDER,
    source_lat=LS_LAT,
    source_lon=LS_LON,
)

if CALCULATE_GF:
    force.compute_greens(modelfile=MODEL_FILE, gfduration=200, T0=-10)
else:
    force.load_greens(modelfile=MODEL_FILE)

force.setup(period_range=PERIOD_RANGE, zeroPhase=True)

#%% INVERT

force.invert(
    zeroTime=119,
    imposeZero=True,
    addtoZero=True,
    jackknife=True,
    num_iter=20,
    frac_delete=0.3,
    alphaset=4.8e-20,
    zeroScaler=2,
    Tikhratio=[0.4, 0.0, 0.6],
)

#%% PLOT

XLIM = (-50, 200)  # [s] x-axis (time) limits for plots
L = 5.8  # [km] Estimate of horizontal COM runout length

# Plot inversion waveform fits and results
force.plotdatafit()
force.plotinv(
    highf_tr=st_hf[0],
    hfshift=hf_shift,
    infra_tr=st_infra[0],
    infra_shift=infra_shift,
    jackshowall=True,
    xlim=XLIM,
    subplots=True,
)
force.plotangmag(xlim=XLIM)

# Calculate/plot trajectories
force.trajectory(
    target_length=L,
    plot_jackknife=True,
    duration=XLIM[1],
    detrend_velocity=XLIM[1]
)
force.trajectory(
    target_length=L,
    plot_jackknife=True,
    duration=XLIM[1],
    detrend_velocity=XLIM[1],
    elevation_profile=True,
)
