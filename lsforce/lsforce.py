import copy
import glob
import os
import pickle
import random as rnd
import subprocess
import tempfile
from shutil import which
from urllib.request import urlretrieve

import numpy as np
import scipy as sp
from matplotlib import lines as mlines
from matplotlib import pyplot as plt
from obspy import Stream, Trace, UTCDateTime, read
from obspy.core import AttribDict
from scipy.signal.windows import triang

# TODO: This is the "beta" URL!
SYNGINE_BASE_URL = 'http://service.iris.washington.edu/iriswsbeta/syngine/1/query?'

# [s] Sampling interval for Green's functions downloaded from Syngine
SYNGINE_DT = 0.25

# [s] Sampling interval for the triangular source time function given to Syngine
TRIANGLE_STF_DT = 0.5

# A nice constant starttime for Syngine and CPS GFs
GF_STARTTIME = UTCDateTime(1900, 1, 1)

# Convert m/s to μm/s
UMS_PER_MS = 1e6


class LSForce:
    r"""Class for performing force inversions.

    Attributes:
        gf_dir (str): Directory containing Green's functions
        gf_computed (bool): Whether or not Green's functions have been computed for this
            object
        filtered_gf_st (:class:`~obspy.core.stream.Stream`): Stream containing filtered
            Green's functions
        inversion_complete (bool): Whether or not the inversion has been run
        filter (dict): Dictionary with keys ``'freqmin'``, ``'freqmax'``,
            ``'zerophase'``, ``'periodmin'``, ``'periodmax'``, and ``'order'``
            specifying filter parameters
        data_length (int): Length in samples of each data trace
        force_sampling_rate (int or float): [Hz] The sampling rate of the force-time
            function
        W (2D array): Weight matrix
        Wvec (1D array): Weight vector
        jackknife (:class:`~obspy.core.util.attribdict.AttribDict`): Dictionary with
            keys ``'Z'``, ``'N'``, ``'E'``, ``'VR_all'``, ``'num_iter'``, and
            ``'frac_delete'`` containing jackknife parameters and results
        angle_magnitude (:class:`~obspy.core.util.attribdict.AttribDict`): Dictionary
            with keys ``'magnitude'``, ``'magnitude_upper'``, ``'magnitude_lower'``,
            ``'vertical_angle'``, and ``'horizontal_angle'`` containing inversion angle
            and magnitude information
        G (2D array): Design matrix
        d (1D array): Data vector
        model: Model vector of concatenated components (n x 1) of solution
        Z: [N] Vertical force time series extracted from model (positive up)
        N: [N] North force time series extracted from model (positive north)
        E: [N] East force time series extracted from model (positive east)
        tvec: [s] Time vector, referenced using `zero_time` (if specified) and corrected
            for `T0` time shift
        VR: [%] Variance reduction. Rule of thumb: This should be ~50–80%, if ~100%,
            solution is fitting data exactly and results are suspect. If ~5%, model may
            be wrong or something else may be wrong with setup
        dtorig: Original data vector
        dtnew: Modeled data vector (Gm-d)
        alpha: Regularization parameter that was used
        alphafit (dict): Dictionary with keys ``'alphas'``, ``'fit'``, and ``'size'``
            specifying regularization parameters tested
    """

    def __init__(
        self, data, data_sampling_rate, nickname=None, main_folder=None, method='full'
    ):
        r"""Create an LSForce object.

        Args:
            data (:class:`~lsforce.lsdata.LSData`): LSData object, corrected for station
                response but not filtered
            data_sampling_rate (int or float): [Hz] Samples per second to use in
                inversion. All data will be resampled to this rate, and Green's
                functions will be created with this rate
            nickname (str): Nickname for this event, used for convenient naming of files
            main_folder (str): If `None`, will use current folder
            method (str): How to parameterize the force-time function. One of `'full'`
                — full waveform inversion using Tikhonov regularization (L2 norm
                minimization) or `'triangle'` — inversion parameterized using
                overlapping triangles, variation on method of Ekström & Stark (2013)
        """

        self.data = data
        self.data_sampling_rate = data_sampling_rate
        self.nickname = nickname
        self.gf_computed = False
        self.inversion_complete = False

        if not main_folder:
            self.main_folder = os.getcwd()
        else:
            self.main_folder = main_folder

        if method not in ['full', 'triangle']:
            raise ValueError(f'Method {method} not yet implemented.')

        self.method = method

    def _get_greens(self):

        # Name the directory where GFs will be stored based on method
        if self.cps_model:
            gf_dir_name = f'cps_{os.path.basename(self.cps_model).split(".")[-2]}'
        else:
            gf_dir_name = f'syngine_{self.syngine_model}'
        self.gf_dir = os.path.join(self.main_folder, gf_dir_name)

        # Label directories containing triangular GFs as such
        if self.method == 'triangle':
            self.gf_dir += f'_triangle_{self.triangle_half_width:g}s'

        # Make GF directory if it doesn't exist
        if not os.path.exists(self.gf_dir):
            os.mkdir(self.gf_dir)

        # Choose the correct delta
        if self.cps_model:
            gf_dt = 1 / self.data_sampling_rate
        else:
            gf_dt = SYNGINE_DT

        # Get list of unique stations
        unique_stations = np.unique([tr.stats.station for tr in self.data.st_proc])

        # Make lists of stations with and without GFs calculated/downloaded
        existing_stations = []
        stations_to_calculate = []
        for station in unique_stations:

            distance = self.data.st_proc.select(station=station)[0].stats.distance
            filename = os.path.join(self.gf_dir, f'{station}.pkl')

            # Check if this EXACT GF exists already
            if os.path.exists(filename):
                stats = read(filename)[0].stats
                if (
                    stats.syngine_model == self.syngine_model
                    and stats.cps_model == self.cps_model
                    and stats.triangle_half_width == self.triangle_half_width
                    and stats.sourcedepthinmeters == self.source_depth
                    and stats.distance == distance
                    and stats.T0 == self.T0
                    and stats.duration == self.gf_duration
                    and stats.delta == gf_dt
                ):
                    gf_exists = True
                else:
                    gf_exists = False
            else:
                gf_exists = False

            # Append to the correct vector
            if gf_exists:
                existing_stations.append(station)
            else:
                stations_to_calculate.append(station)

        # Initalize empty Stream to hold all GFs
        st_gf = Stream()

        # CPS
        if self.cps_model:

            # If we have to calculate some stations, go through the process
            if stations_to_calculate:

                # Print status
                print(
                    f'Calculating Green\'s functions for {len(stations_to_calculate)} '
                    'station(s):'
                )
                for station in stations_to_calculate:
                    print(f'\t{station}')

                # Create a temporary directory to run CPS in, and change into it
                cwd = os.getcwd()
                temp_dir = tempfile.TemporaryDirectory()
                os.chdir(temp_dir.name)

                # Write the "dist" file
                gf_length_samples = int(self.gf_duration * self.data_sampling_rate)
                with open('dist', 'w') as f:
                    for sta in stations_to_calculate:
                        dist = self.data.st_proc.select(station=sta)[0].stats.distance
                        f.write(f'{dist} {gf_dt} {gf_length_samples} {self.T0} 0\n')

                # Run hprep96 and hspec96
                subprocess.call(
                    [
                        'hprep96',
                        '-HR',
                        '0.',
                        '-HS',
                        str(self.source_depth / 1000),  # Converting to km for CPS
                        '-M',
                        self.cps_model,
                        '-d',
                        'dist',
                        '-R',
                        '-EXF',
                    ]
                )
                with open('hspec96.out', 'w') as f:
                    subprocess.call('hspec96', stdout=f)

                # Run hpulse96 (using the multiplier here to get a 1 N impulse), also
                # keep track of pulse half-width so we can make the GFs acausal later
                args = ['hpulse96', '-d', 'dist', '-m', f'{1e-10:.10f}', '-OD', '-V']
                if self.triangle_half_width is not None:
                    args += ['-t', '-l', str(int(self.triangle_half_width / gf_dt))]
                    pulse_half_width = self.triangle_half_width  # [s]
                else:
                    args += ['-p']
                    pulse_half_width = 2 * gf_dt  # [s]
                with open('Green', 'w') as f:
                    subprocess.call(args, stdout=f)

                # Convert to SAC files
                subprocess.call(['f96tosac', 'Green'])

                # Go through and read in files (same order as dist file)
                for i, station in enumerate(stations_to_calculate):
                    for file in glob.glob(f'B{i + 1:03d}1???F.sac'):

                        # Grab stats of data trace
                        stats = self.data.st_proc.select(station=station)[0].stats

                        # Read GF in as an ObsPy Trace
                        gf_tr = read(file)[0]

                        # Add metadata
                        gf_tr.stats.network = stats.network
                        gf_tr.stats.station = station
                        gf_tr.stats.location = 'SE'
                        gf_tr.stats.distance = stats.distance
                        gf_tr.stats.syngine_model = self.syngine_model
                        gf_tr.stats.cps_model = self.cps_model
                        gf_tr.stats.sourcedepthinmeters = self.source_depth
                        gf_tr.stats.T0 = self.T0
                        gf_tr.stats.duration = self.gf_duration
                        gf_tr.stats.triangle_half_width = self.triangle_half_width

                        gf_tr.data *= 0.01  # Convert from cm to m

                        # Add Trace to overall GF Stream
                        st_gf += gf_tr

                # Clean up
                temp_dir.cleanup()
                os.chdir(cwd)

                # Trim to length (gf_duration - T0) seconds, and correct for the pulse
                # half-width to make GFs acausal (since the -Z flag doesn't work!)
                starttime = st_gf[0].stats.starttime
                endtime = starttime + self.gf_duration - self.T0
                st_gf.trim(starttime + pulse_half_width, endtime + pulse_half_width)

                # Give a nice start time
                for tr in st_gf:
                    tr.stats.starttime = GF_STARTTIME

                # Save as individual files
                for station in stations_to_calculate:
                    filename = os.path.join(self.gf_dir, f'{station}.pkl')
                    st_gf.select(station=station).write(filename, format='PICKLE')

            # Now just load in the GFs which already exist
            for i, station in enumerate(existing_stations):
                filename = os.path.join(self.gf_dir, f'{station}.pkl')
                st_gf += read(filename)

                # Print status
                progress = len(stations_to_calculate) + i + 1
                print(f'Found {station} ({progress}/{len(unique_stations)})')

        # Syngine
        else:

            # Go station-by-station
            for i, station in enumerate(unique_stations):

                filename = os.path.join(self.gf_dir, f'{station}.pkl')

                # Either load this GF if it already exists, or download it
                if station in existing_stations:
                    st_syn = read(filename)
                else:
                    # Grab stats of data trace
                    stats = self.data.st_proc.select(station=station)[0].stats

                    # Get GFs for this station
                    st_syn = self._get_greens_for_station(
                        network=stats.network,
                        station=station,
                        back_azimuth=stats.back_azimuth,
                        distance=stats.distance,
                    )

                    # TODO: Understand why we have to invert these!
                    for channel in 'RHF', 'RVF', 'THF':
                        for tr in st_syn.select(channel=channel):
                            tr.data *= -1

                    st_syn.write(filename, format='PICKLE')

                # Add this station's GF's to overall Stream
                st_gf += st_syn

                # Print status
                if station in existing_stations:
                    action_string = 'Found'
                else:
                    action_string = 'Downloaded'
                print(f'{action_string} {station} ({i + 1}/{len(unique_stations)})')

        self.gf_computed = True

        return st_gf

    def _get_greens_for_station(self, network, station, back_azimuth, distance):

        # Provide triangle STF params if we're using the triangle method
        if self.method == 'triangle':
            stf_offset = self.triangle_half_width  # Ensure peak of triangle at t=0
            stf_spacing = TRIANGLE_STF_DT
            # The below construction ensures the triangle is centered on 1 and goes to 0
            # at each end, e.g. [0, 0.5, 1, 0.5, 0] instead of [0.25, 0.75, 0.75, 0.25]
            stf_data = np.hstack(
                [
                    0,
                    triang((int(self.triangle_half_width / TRIANGLE_STF_DT) * 2) - 1),
                    0,
                ]
            )
            read_func = _read  # Use the long-URL wrapper for ObsPy read
        else:
            stf_offset = None
            stf_spacing = None
            stf_data = None
            read_func = read  # Just directly read using ObsPy

        # Convert to radians for NumPy
        back_azimuth_radians = np.deg2rad(back_azimuth)

        # Vertical force (downward)
        st_vf = read_func(
            self._build_syngine_url(
                network=network,
                station=station,
                components='ZR',
                forces=(-1, 0, 0),
                stf_offset=stf_offset,
                stf_spacing=stf_spacing,
                stf_data=stf_data,
            )
        )
        for tr in st_vf.select(component='Z'):
            tr.stats.channel = 'ZVF'
        for tr in st_vf.select(component='R'):
            tr.stats.channel = 'RVF'

        # Horizontal force (radial)
        st_hf_r = read_func(
            self._build_syngine_url(
                network=network,
                station=station,
                components='ZR',
                forces=(0, np.cos(back_azimuth_radians), -np.sin(back_azimuth_radians)),
                stf_offset=stf_offset,
                stf_spacing=stf_spacing,
                stf_data=stf_data,
            )
        )
        for tr in st_hf_r.select(component='Z'):
            tr.stats.channel = 'ZHF'
        for tr in st_hf_r.select(component='R'):
            tr.stats.channel = 'RHF'

        # Horizontal force (transverse)
        st_hf_t = read_func(
            self._build_syngine_url(
                network=network,
                station=station,
                components='T',
                forces=(
                    0,
                    -np.sin(back_azimuth_radians),
                    -np.cos(back_azimuth_radians),
                ),
                stf_offset=stf_offset,
                stf_spacing=stf_spacing,
                stf_data=stf_data,
            )
        )
        for tr in st_hf_t.select(component='T'):
            tr.stats.channel = 'THF'

        # Assemble big Stream
        st_syn = st_vf + st_hf_r + st_hf_t

        # Add metadata and sort
        for tr in st_syn:
            tr.stats.distance = distance
            tr.stats.cps_model = self.cps_model
            tr.stats.syngine_model = self.syngine_model
            tr.stats.sourcedepthinmeters = self.source_depth
            tr.stats.T0 = self.T0
            tr.stats.duration = self.gf_duration
            tr.stats.triangle_half_width = self.triangle_half_width
            tr.stats.starttime = GF_STARTTIME  # Give a nice start time
        st_syn.sort(keys=['channel'])

        return st_syn

    def _build_syngine_url(
        self,
        network,
        station,
        components,
        forces,
        stf_offset=None,
        stf_spacing=None,
        stf_data=None,
    ):

        parameters = [
            'format=miniseed',
            'components=' + components,
            'units=displacement',
            'model=' + self.syngine_model,
            'dt=' + str(SYNGINE_DT),
            'starttime=' + str(self.T0),
            'endtime=' + str(self.gf_duration - self.T0),
            'network=' + network,
            'station=' + station,
            'sourcelatitude=' + str(self.data.source_lat),
            'sourcelongitude=' + str(self.data.source_lon),
            'sourcedepthinmeters=' + str(self.source_depth),
            'sourceforce=' + ','.join([str(force) for force in forces]),
            'nodata=404',
        ]

        if stf_offset is not None and stf_spacing is not None and stf_data is not None:
            parameters += [
                'cstf-relative-origin-time-in-sec=' + str(stf_offset),
                'cstf-sample-spacing-in-sec=' + str(stf_spacing),
                'cstf-data=' + ','.join([str(sample) for sample in stf_data]),
            ]
        elif stf_offset is not None or stf_spacing is not None or stf_data is not None:
            raise ValueError('All three CSTF parameters must be provided!')

        url = SYNGINE_BASE_URL + '&'.join(parameters)

        return url

    def setup(
        self,
        period_range,
        syngine_model=None,
        cps_model=None,
        triangle_half_width=None,
        source_depth=0,
        weights=None,
        noise_window_dur=None,
        filter_order=2,
        zerophase=False,
    ):
        r"""Downloads/computes Green's functions (GFs) and creates all matrices.

        Args:
            period_range (list or tuple): [s] Bandpass filter corners
            syngine_model (str): Name of Syngine model to use. If this is not None, then
                we calculate GFs using Syngine (preferred)
            cps_model (str): Filename of CPS model to use. If this is not None, then we
                calculate GFs using CPS
            triangle_half_width (int or float): [s] Half-width of triangles; only used
                if the triangle method is being used
            source_depth (int or float): [m] Source depth in meters
            weights (list or tuple or str): If `None`, no weighting is applied. An array
                of floats with length ``st.count()`` and in the order of the `st`
                applies manual weighting. If `'prenoise'`, uses standard deviation of
                a noise window defined by `noise_window_dur` to weight. If `'distance'`,
                weights by 1/distance
            noise_window_dur (int or float): [s] Length of noise window for `'prenoise'`
                weighting scheme (if not `None`, `weights` is set to `'prenoise'`)
            filter_order (int): Order of filter applied over period_range
            zerophase (bool): If `True`, zero-phase filtering will be used
        """

        self.syngine_model = syngine_model
        self.source_depth = source_depth

        # Explicitly ignore the triangle half-width parameter if it's not relevant
        if self.method != 'triangle' and triangle_half_width is not None:
            triangle_half_width = None
            print(
                'Ignoring `triangle_half_width` parameter since you\'re not using the '
                'triangle method.'
            )

        # Make sure user specifies the triangle half-width if they want that method
        if self.method == 'triangle' and triangle_half_width is None:
            raise ValueError('triangle method is specified but no half-width given!')
        self.triangle_half_width = triangle_half_width

        # If user wants CPS to be run, make sure that 1) they have it installed and 2)
        # they have provided a valid filepath
        if cps_model:
            # 1) Is CPS installed?
            if not which('hprep96'):
                raise OSError(
                    'CPS Green\'s function calculation requested, but CPS not found on '
                    'system. Install CPS and try again.'
                )
            # 2) Is `cps_model` a file?
            if not os.path.exists(cps_model):
                raise OSError(f'Could not find CPS model file "{cps_model}"')
            else:
                cps_model = os.path.abspath(cps_model)  # Get full path
        self.cps_model = cps_model

        # The user must specify ONE of `syngine_model` and `cps_model`
        if (self.syngine_model and self.cps_model) or (
            not self.syngine_model and not self.cps_model
        ):
            raise ValueError('You must specify ONE of syngine_model or cps_model!')

        # Automatically choose an appropriate T0 and GF duration based on data/method
        if self.method == 'triangle':
            self.T0 = -2 * self.triangle_half_width  # [s] Double the half-width
        else:
            self.T0 = -10  # [s]
        min_time = np.min([tr.stats.starttime for tr in self.data.st_proc])
        max_time = np.max([tr.stats.endtime for tr in self.data.st_proc])
        self.gf_duration = max_time - min_time  # [s]

        # Create filter dictionary to keep track of filter used without creating too
        # many new attributes
        self.filter = {
            'freqmin': 1.0 / period_range[1],
            'freqmax': 1.0 / period_range[0],
            'zerophase': zerophase,
            'periodmin': period_range[0],
            'periodmax': period_range[1],
            'order': filter_order,
        }

        # Clear weights
        self.Wvec = None
        self.W = None

        if weights is None:
            # Don't weight at all
            weight_method = None
        elif isinstance(weights, str):
            # The user specified a weight method
            weight_method = weights
        else:
            # The user specified a vector of station weights
            weight_method = 'Manual'
            self.weights = weights

        if weight_method != 'Manual':
            if weights == 'prenoise' and noise_window_dur is None:
                raise ValueError(
                    'noise_window_dur must be defined if prenoise weighting is used.'
                )

        # check if sampling rate specified is compatible with period_range
        if 2.0 * self.filter['freqmax'] > self.data_sampling_rate:
            raise ValueError(
                'data_sampling_rate and period_range are not compatible (violates '
                'Nyquist).'
            )

        # Always work on copy of data
        st = self.data.st_proc.copy()

        # filter data to band specified
        st.filter(
            'bandpass',
            freqmin=self.filter['freqmin'],
            freqmax=self.filter['freqmax'],
            corners=self.filter['order'],
            zerophase=self.filter['zerophase'],
        )

        # resample st to data_sampling_rate
        st.resample(self.data_sampling_rate)

        # make sure st data are all the same length
        lens = [len(trace.data) for trace in st]
        if len(set(lens)) != 1:
            print(
                'Resampled records are of differing lengths. Interpolating all records '
                'to same start time and sampling rate.'
            )
            stts = [tr.stats.starttime for tr in st]
            lens = [tr.stats.npts for tr in st]
            st.interpolate(
                self.data_sampling_rate, starttime=np.max(stts), npts=np.min(lens) - 1
            )

        self.data_length = st[0].stats.npts
        total_data_length = self.data_length * st.count()

        # Load in GFs
        print('Getting Green\'s functions...')
        st_gf = self._get_greens()

        # Process GFs in bulk
        st_gf.detrend()
        st_gf.taper(max_percentage=0.05)
        st_gf.filter(
            'bandpass',
            freqmin=self.filter['freqmin'],
            freqmax=self.filter['freqmax'],
            corners=self.filter['order'],
            zerophase=self.filter['zerophase'],
        )
        if self.syngine_model:  # Only need to do this if Syngine
            st_gf.interpolate(
                sampling_rate=self.data_sampling_rate, method='lanczos', a=20
            )
        st_gf.sort(keys=['channel'])

        # Store the filtered GFs
        self.filtered_gf_st = st_gf

        # Initialize weighting matrices
        Wvec = np.ones(total_data_length)
        indx = 0
        weight = np.ones(self.data.st_proc.count())

        # Store data length
        n = self.data_length

        if self.method == 'full':
            # Set sampling rate
            self.force_sampling_rate = self.data_sampling_rate
        elif self.method == 'triangle':
            self.force_sampling_rate = 1.0 / self.triangle_half_width
            # Number of samples to shift each triangle by
            fshiftby = int(self.triangle_half_width * self.data_sampling_rate)
            # Number of shifts, corresponds to length of force time function
            Flen = int(np.floor(self.data_length / fshiftby))
            # Triangle GFs are multiplied by the triangle half-width so that they
            # reflect the ground motion induced for a triangle with PEAK 1 N instead of
            # AREA of 1 N*s
            for tr in st_gf:
                tr.data = tr.data * self.triangle_half_width
        else:
            raise ValueError(f'Method {self.method} not supported.')

        for i, tr in enumerate(st):

            # Find component and station of Trace
            component = tr.stats.channel[-1]
            station = tr.stats.station

            if component == 'Z':
                zvf = st_gf.select(station=station, channel='ZVF')[0]
                zhf = st_gf.select(station=station, channel='ZHF')[0]
                if self.method == 'full':
                    ZVF = _makeconvmat(zvf.data, size=(n, n))
                    ZHF = _makeconvmat(zhf.data, size=(n, n))
                else:
                    ZVF = _makeshiftmat(zvf.data, shiftby=fshiftby, size1=(n, Flen))
                    ZHF = _makeshiftmat(zhf.data, shiftby=fshiftby, size1=(n, Flen))
                az_radians = np.deg2rad(tr.stats.azimuth)
                newline = np.hstack(
                    [ZVF, ZHF * np.cos(az_radians), ZHF * np.sin(az_radians)]
                )

            elif component == 'R':
                rvf = st_gf.select(station=station, channel='RVF')[0]
                rhf = st_gf.select(station=station, channel='RHF')[0]
                if self.method == 'full':
                    RVF = _makeconvmat(rvf.data, size=(n, n))
                    RHF = _makeconvmat(rhf.data, size=(n, n))
                else:
                    RVF = _makeshiftmat(rvf.data, shiftby=fshiftby, size1=(n, Flen))
                    RHF = _makeshiftmat(rhf.data, shiftby=fshiftby, size1=(n, Flen))
                az_radians = np.deg2rad(tr.stats.azimuth)
                newline = np.hstack(
                    [RVF, RHF * np.cos(az_radians), RHF * np.sin(az_radians)]
                )

            elif component == 'T':
                thf = st_gf.select(station=station, channel='THF')[0]
                if self.method == 'full':
                    THF = _makeconvmat(thf.data, size=(n, n))
                else:
                    THF = _makeshiftmat(thf.data, shiftby=fshiftby, size1=(n, Flen))
                TVF = 0.0 * THF.copy()  # Just zeros for TVF
                az_radians = np.deg2rad(tr.stats.azimuth)
                newline = np.hstack(
                    [TVF, THF * np.sin(az_radians), -THF * np.cos(az_radians)]
                )

            else:
                raise ValueError(f'Data not rotated to ZRT for {station}.')

            # Deal with data
            datline = tr.data

            if i == 0:  # initialize G and d if first station
                G = newline.copy()
                d = datline.copy()
            else:  # otherwise build on G and d
                G = np.vstack((G, newline.copy()))
                d = np.hstack((d, datline.copy()))

            if weights is not None:
                if weight_method == 'Manual':
                    weight[i] = weights[i]
                elif weights == 'prenoise':
                    weight[i] = 1.0 / np.std(
                        tr.data[0 : int(noise_window_dur * tr.stats.sampling_rate)]
                    )
                elif weights == 'distance':
                    weight[i] = tr.stats.distance

                Wvec[indx : indx + self.data_length] = (
                    Wvec[indx : indx + self.data_length] * weight[i]
                )
                indx += self.data_length
        if self.method == 'full':
            # Need to multiply G by sample interval [s] since convolution is an integral
            self.G = G * 1.0 / self.data_sampling_rate
        else:
            # We don't need to scale the triangle method GFs by the sample rate since
            # this method is not a convolution
            self.G = G

        # Normalize Wvec so largest weight is 1.
        self.Wvec = Wvec / np.max(np.abs(Wvec))
        self.weights = weight / np.max(np.abs(weight))

        if np.shape(G)[0] != len(d):
            raise ValueError(
                'G and d sizes are not compatible, fix something somewhere.'
            )

        self.d = d
        if weights is not None:
            self.W = np.diag(self.Wvec)
        else:
            self.W = None

    def invert(
        self,
        zero_time=None,
        impose_zero=False,
        add_to_zero=False,
        max_duration=None,
        jackknife=False,
        num_iter=200,
        frac_delete=0.5,
        **kwargs,
    ):
        r"""Performs single-force inversion using Tikhonov regularization.

        Args:
            zero_time (int or float): [s] Optional estimated start time of real
                (avalanche-related) part of signal, in seconds from start time of
                seismic data. Useful for making figures showing selected start time and
                also for the `impose_zero` option
            impose_zero (bool): Adds weighting matrix to suggest that forces tend
                towards zero prior to `zero_time` (`zero_time` must be defined)
            add_to_zero (bool): Adds weighting matrix to suggest that all components of
                force integrate to zero
            max_duration (int or float): Maximum duration allowed for the event,
                starting at `zero_time` if defined, otherwise starting from the
                beginning of the seismic data. Forces after this will tend towards zero.
                This helps tamp down artifacts due to edge effects, etc.
            jackknife (bool): If `True`, perform `num_iter` additional iterations of the
                model while randomly discarding `frac_delete` of the data
            num_iter (int): Number of jackknife iterations to perform
            frac_delete (int or float): Fraction (out of 1) of data to discard for each
                iteration
            **kwargs: Additional keyword arguments to be passed on to the inversion
                method
        """

        # Check inputs for consistency
        if impose_zero and not zero_time:
            raise ValueError('impose_zero set to True but no zero_time provided.')

        # Save input choices
        self.add_to_zero = add_to_zero
        self.zero_time = zero_time
        self.impose_zero = impose_zero
        self.max_duration = max_duration

        # Initialize (also serves to clear any previous results if this is a rerun)
        self.model = None
        self.Z = None
        self.N = None
        self.E = None
        self.tvec = None
        self.VR = None  # variance reduction
        self.dtorig = None  #
        self.dtnew = None
        self.alpha = None
        self.alphafit = {'alphas': None, 'fit': None, 'size': None}

        if jackknife:
            # Initialize
            self.jackknife = AttribDict(
                Z=AttribDict(lower=[], upper=[], all=[]),
                N=AttribDict(lower=[], upper=[], all=[]),
                E=AttribDict(lower=[], upper=[], all=[]),
                VR_all=[],
                num_iter=num_iter,
                frac_delete=frac_delete,
            )
        else:
            self.jackknife = None

        self._tikinvert(**kwargs)

        self.inversion_complete = True

    def _tikinvert(
        self,
        alphaset=None,
        zero_scaler=2.0,
        zero_taper_length=20.0,
        tikhonov_ratios=(1.0, 0.0, 0.0),
    ):
        r"""Performs full-waveform inversion using Tikhonov regularization.

        Args:
            alphaset (int or float): Set regularization parameter. If `None`, will
                search for best alpha using the L-curve method
            zero_scaler (int or float): Relative strength of zero constraint from
                0 to 10. The lower the number, the weaker
                the constraint. Values up to 30 are technically allowed but
                discouraged because high zero_scaler values risk
                the addition of high frequency oscillations due to the sudden release
                of the constraint.
            zero_taper_length (int or float): [s] Length of taper for `zero_scaler`.
                Tapers that are toos hort can result in sharp spiky artifacts.
            tikhonov_ratios (list or tuple): Proportion each regularization method
                contributes to the overall regularization effect, where values
                correspond to [0th order, 1st order, 2nd order]. Must sum to 1
        """
        if zero_scaler < 0.0 or zero_scaler > 30.0:
            raise ValueError('zero_scaler cannot be less than 0 or more than 30')

        if np.sum(tikhonov_ratios) != 1.0:
            raise ValueError('Tikhonov ratios must add to 1.')

        if self.W is not None:
            Ghat = self.W.dot(self.G)
            dhat = self.W.dot(self.d)
        else:
            Ghat = self.G
            dhat = self.d

        if self.jackknife is not None:  # save version at this point for use later
            Ghatori = Ghat.copy()
            dhatori = dhat.copy()

        m, n = np.shape(Ghat)

        Ghatnorm = np.linalg.norm(Ghat)

        dl = self.data_length
        gl = int(n / 3)  # Force vector length

        if self.add_to_zero is True:  # constrain forces to add to zero
            scaler = Ghatnorm
            first1 = np.hstack((np.ones(gl), np.zeros(2 * gl)))
            second1 = np.hstack((np.zeros(gl), np.ones(gl), np.zeros(gl)))
            third1 = np.hstack((np.zeros(2 * gl), np.ones(gl)))
            A1 = np.vstack((first1, second1, third1)) * scaler
            Ghat = np.vstack((Ghat, A1))
            dhat = np.hstack((dhat, np.zeros(3)))
        else:
            A1 = None

        scaler = Ghatnorm * (zero_scaler / 30.0)

        if self.impose_zero:  # tell model when there should be no forces
            len2 = int(((self.zero_time + self.T0) * self.force_sampling_rate))
            if self.method == 'triangle':
                # No taper
                vals2 = np.hstack((np.ones(len2), np.zeros(gl - len2)))
            elif self.method == 'full':
                # Taper
                len3 = int(
                    zero_taper_length * self.force_sampling_rate
                )  # make it constant
                temp = np.hanning(2 * len3)
                temp = temp[len3:]
                vals2 = np.hstack((np.ones(len2 - len3), temp))
            else:
                raise NotImplementedError

            for i, val in enumerate(vals2):
                first1 = np.zeros(3 * gl)
                second1 = first1.copy()
                third1 = first1.copy()
                if i > 0:
                    first1[i] = val
                    second1[i + gl] = val
                    third1[i + 2 * gl] = val
                if i == 0:
                    A2 = np.vstack((first1, second1, third1))
                else:
                    A2 = np.vstack((A2, first1, second1, third1))
            A2 *= scaler
            Ghat = np.vstack((Ghat, A2))
            dhat = np.hstack((dhat, np.zeros(len(vals2) * 3)))
        else:
            A2 = None

        if self.max_duration is not None:
            if self.zero_time is None:
                zerotime = 0.0
            else:
                zerotime = self.zero_time
            startind = int(
                (zerotime + self.T0 + self.max_duration) * self.force_sampling_rate
            )
            if self.method == 'triangle':
                vals3 = np.zeros(gl)
                vals3[startind:] = 1.0  # no taper
                for i, val in enumerate(vals3):
                    first1 = np.zeros(3 * gl)
                    second1 = first1.copy()
                    third1 = first1.copy()
                    if val > 0:
                        first1[i] = val
                        second1[i + gl] = val
                        third1[i + 2 * gl] = val
                    if i == 0:
                        A3 = np.vstack((first1, second1, third1))
                    else:
                        A3 = np.vstack((A3, first1, second1, third1))
            if self.method == 'full':
                len2 = int(gl - startind)
                len3 = int(
                    np.round(0.2 * len2)
                )  # 20% taper so zero imposition isn't sudden
                temp = np.hanning(2 * len3)
                temp = temp[:len3]
                vals3 = np.hstack((temp, np.ones(len2 - len3)))

                for i, val in enumerate(vals3):
                    place = i + startind
                    first1 = np.zeros(3 * gl)
                    second1 = first1.copy()
                    third1 = first1.copy()
                    if val > 0:
                        first1[place] = val
                        second1[place + gl] = val
                        third1[place + 2 * gl] = val
                    if i == 0:
                        A3 = np.vstack((first1, second1, third1))
                    else:
                        A3 = np.vstack((A3, first1, second1, third1))
            A3 *= scaler

            Ghat = np.vstack((Ghat, A3))
            dhat = np.hstack((dhat, np.zeros(len(vals3) * 3)))
        else:
            A3 = None

        if alphaset is not None:
            alpha = alphaset
        dhat = dhat.T

        # Build roughening matrix
        I = np.eye(n, n)
        if tikhonov_ratios[1] != 0.0:
            # Build L1 (first order) roughening matrix
            L1 = np.diag(-1 * np.ones(n)) + np.diag(np.ones(n - 1), k=1)
            L1part = L1.T @ L1
        else:
            L1part = 0.0
            L1 = 0.0
        if tikhonov_ratios[2] != 0.0:
            # Build L2 (second order) roughening matrix
            L2 = (
                np.diag(np.ones(n))
                + np.diag(-2 * np.ones(n - 1), k=1)
                + np.diag(np.ones(n - 2), k=2)
            )
            L2part = L2.T @ L2
        else:
            L2 = 0.0
            L2part = 0.0

        if alphaset is None:
            alpha, fit1, size1, alphas = _find_alpha(
                Ghat, dhat, I, L1, L2, tikhonov_ratios=tikhonov_ratios
            )
            print(f'best alpha is {alpha:6.1e}')
            self.alpha = alpha
            self.alphafit['alphas'] = alphas
            self.alphafit['fit'] = fit1
            self.alphafit['size'] = size1
        else:
            self.alpha = alpha

        Apart = Ghat.conj().T @ Ghat

        A = Apart + alpha ** 2 * (
            tikhonov_ratios[0] * I
            + tikhonov_ratios[1] * L1part
            + tikhonov_ratios[2] * L2part
        )  # Combo of all regularization things (if any are zero they won't matter)
        x = np.squeeze(Ghat.conj().T @ dhat)

        model, residuals, rank, s = sp.linalg.lstsq(A, x)
        self.model = model.copy()
        div = int(len(model) / 3)

        # Flip so up is positive
        self.Z = -model[0:div]
        self.N = model[div : 2 * div]
        self.E = model[2 * div :]

        dtnew = self.G.dot(model)
        self.dtnew = np.reshape(dtnew, (self.data.st_proc.count(), dl))
        self.dtorig = np.reshape(self.d, (self.data.st_proc.count(), dl))

        # Compute variance reduction
        # TODO compute only for where zeroing not applied (use taper)
        self.VR = _varred(self.dtorig, self.dtnew)
        print(f'Variance reduction = {self.VR:f} percent')
        tvec = (
            np.arange(
                0,
                len(self.Z) * 1 / self.force_sampling_rate,
                1 / self.force_sampling_rate,
            )
            - self.T0
        )

        # Adjust time vector for zero time, if one was provided
        if self.zero_time is not None:
            tvec -= self.zero_time

        self.tvec = tvec
        self.dtvec = np.arange(
            0, self.data_length / self.data_sampling_rate, 1 / self.data_sampling_rate
        )
        if self.zero_time is not None:
            self.dtvec -= self.zero_time
        # Use constant alpha parameter (found above, if not previously set) for
        # jackknife iterations
        stasets = []
        if self.jackknife is not None:
            # Start jackknife iterations
            for ii in range(self.jackknife.num_iter):
                numcut = int(
                    round(self.jackknife.frac_delete * self.data.st_proc.count())
                )
                numkeep = self.data.st_proc.count() - numcut
                indxcut = rnd.sample(list(range(self.data.st_proc.count())), numcut)
                stasets.append(indxcut)

                obj = [
                    sum(ind)
                    for ind in zip(
                        np.tile(list(range(self.data_length)), len(indxcut)),
                        np.repeat(
                            [x1 * self.data_length for x1 in indxcut], self.data_length
                        ),
                    )
                ]

                dhat1 = np.delete(dhatori.copy(), obj)
                Ghat1 = np.delete(Ghatori.copy(), obj, axis=0)

                Gtemp = np.delete(self.G.copy(), obj, axis=0)
                dtemp = np.delete(self.d.copy(), obj)

                if A1 is not None:  # Use A1, A2, A3 from full solution, if exist
                    Ghat1 = np.vstack((Ghat1, A1))
                    dhat1 = np.hstack((dhat1, np.zeros(3)))
                if A2 is not None:
                    Ghat1 = np.vstack((Ghat1, A2))
                    dhat1 = np.hstack((dhat1, np.zeros(len(vals2) * 3)))
                if A3 is not None:
                    Ghat1 = np.vstack((Ghat1, A3))
                    dhat1 = np.hstack((dhat1, np.zeros(len(vals3) * 3)))

                dhat1 = dhat1.T
                Apart = Ghat1.conj().T @ Ghat1

                # Combo of all regularization things (if any are zero they won't matter)
                Aj = Apart + self.alpha ** 2 * (
                    tikhonov_ratios[0] * I
                    + tikhonov_ratios[1] * L1part
                    + tikhonov_ratios[2] * L2part
                )
                xj = np.squeeze(Ghat1.conj().T @ dhat1)

                model, residuals, rank, s = sp.linalg.lstsq(Aj, xj)
                div = int(len(model) / 3)
                Zf = -model[0:div]  # Flip so up is positive
                Nf = model[div : 2 * div]
                Ef = model[2 * div :]
                dtnew = Gtemp.dot(model)
                dt = np.reshape(dtemp, (numkeep, dl))

                VR = _varred(dt, dtnew)
                self.jackknife.Z.all.append(Zf.copy())
                self.jackknife.N.all.append(Nf.copy())
                self.jackknife.E.all.append(Ef.copy())
                self.jackknife.VR_all.append(VR.copy())

            self.jackknife.Z.lower = np.percentile(self.jackknife.Z.all, 2.5, axis=0)
            self.jackknife.Z.upper = np.percentile(self.jackknife.Z.all, 97.5, axis=0)
            self.jackknife.E.lower = np.percentile(self.jackknife.E.all, 2.5, axis=0)
            self.jackknife.E.upper = np.percentile(self.jackknife.E.all, 97.5, axis=0)
            self.jackknife.N.lower = np.percentile(self.jackknife.N.all, 2.5, axis=0)
            self.jackknife.N.upper = np.percentile(self.jackknife.N.all, 97.5, axis=0)

            self.jackknife.VR_all = np.array(self.jackknife.VR_all)

            print(
                f'Jackknife VR stats: max {self.jackknife.VR_all.max():2.0f}, '
                f'min {self.jackknife.VR_all.min():2.0f}, '
                f'median {np.median(self.jackknife.VR_all):2.0f}'
            )

    def plot_fits(self, equal_scale=True, xlim=None):
        r"""Create a plot showing the model-produced waveform fit to the data.

        Args:
            equal_scale (bool): If `True`, all plots will share the same y-axis scale
            xlim (list or tuple): [s] Array (length two) of x-axis limits (time relative
                to zero time)

        Returns:
            :class:`~matplotlib.figure.Figure`: Output figure handle
        """

        data_color = 'black'
        model_color = 'red'
        lw = 1

        if not equal_scale:
            amp_string = 'normalized'
            spacing = 2
        else:
            amp_string = 'equal scale'
            spacing = np.abs(self.dtorig).max() * 2  # Spacing controlled by largest amp

        fig, ax = plt.subplots(figsize=(8, 12))

        offset = 0
        yticks = []
        yticklabels = []

        for i, tr in enumerate(self.data.st_proc):
            if not equal_scale:
                scaler = np.abs(self.dtorig[i].max())  # Scale to original data
            else:
                scaler = 1
            ax.plot(
                self.dtvec,
                self.dtorig[i] / scaler + offset,
                color=data_color,
                linewidth=lw,
            )
            ax.plot(
                self.dtvec,
                self.dtnew[i] / scaler + offset,
                color=model_color,
                linewidth=lw,
            )
            label = (
                f'{tr.stats.network}.{tr.stats.station} ({tr.stats.channel[-1]}) '
                f'– {tr.stats.distance:.1f} km'
            )
            yticklabels.append(label)
            yticks.append(offset)
            offset -= spacing

        # Misc. tweaks
        if xlim:
            ax.set_xlim(xlim)
        else:
            ax.set_xlim(self.dtvec[0], self.dtvec[-1])
        ax.set_ylim(yticks[-1] - spacing / 2, yticks[0] + spacing / 2)
        ax.set_xlabel('Time (s)')
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticklabels)
        ax.tick_params(left=False, top=True, which='both')
        ax.yaxis.set_tick_params(length=0, pad=4)

        ax.set_title(f'Variance reduction {self.VR:.1f}% ({amp_string})', pad=10)

        # Make legend
        data_handle = mlines.Line2D(
            [], [], color=data_color, label='Data', linewidth=lw
        )
        model_handle = mlines.Line2D(
            [], [], color=model_color, label='Model', linewidth=lw
        )
        ax.legend(handles=[data_handle, model_handle], loc='upper right')

        fig.tight_layout()
        fig.show()

        return fig

    def plot_forces(
        self,
        subplots=False,
        xlim=None,
        ylim=None,
        sameY=True,
        highf_tr=None,
        hfshift=0.0,
        hfylabel=None,
        infra_tr=None,
        infra_shift=0,
        tvecshift=0.0,
        jackshowall=False,
    ):
        r"""Plot inversion result.

        Args:
            subplots (bool): If `True`, make subplots for components, otherwise plot all
                on one plot
            xlim (list or tuple): x-axis limits
            ylim (list or tuple): y-axis limits
            sameY (bool): If `True`, use same y-axis limits for all plots
            highf_tr (:class:`~obspy.core.trace.Trace`): Seismic trace with start time
                identical to start time of the data used in the inversion
            hfshift (int or float): [s] Time shift for seismic trace
            hfylabel (str): Label used for seismic trace. If not defined, will use
                station name
            infra_tr (:class:`~obspy.core.trace.Trace`): Infrasound trace with start
                time identical to start time of the data used in the inversion
            infra_shift (int or float): [s] Time shift for infrasound trace
            tvecshift (int or float): [s] Shift time vector manually by this many
                seconds, not usually needed, only for display purposes (REMOVE?)
            jackshowall (bool): If `True` and jackknife was run, will show all
                individual runs (changes `subplots` to `True`)

        Returns:
            :class:`~matplotlib.figure.Figure`: Output figure handle
        """

        tvec = self.tvec - tvecshift

        annot_kwargs = dict(xy=(0.99, 0.25), xycoords='axes fraction', ha='right')

        # Find y limits
        if self.jackknife is None:
            if ylim is None:
                ylim1 = (
                    np.amin([self.Z.min(), self.E.min(), self.N.min()]),
                    np.amax([self.Z.max(), self.E.max(), self.N.max()]),
                )
                ylim = (
                    ylim1[0] + 0.1 * ylim1[0],
                    ylim1[1] + 0.1 * ylim1[1],
                )  # add 10% on each side to make it look nicer
        else:
            Zupper = self.jackknife.Z.upper
            Nupper = self.jackknife.N.upper
            Eupper = self.jackknife.E.upper
            Zlower = self.jackknife.Z.lower
            Nlower = self.jackknife.N.lower
            Elower = self.jackknife.E.lower
            if ylim is None:
                ylim1 = (
                    np.amin(
                        [
                            Zlower.min(),
                            Elower.min(),
                            Nlower.min(),
                            self.Z.min(),
                            self.N.min(),
                            self.E.min(),
                        ]
                    ),
                    np.amax(
                        [
                            Zupper.max(),
                            Eupper.max(),
                            Nupper.max(),
                            self.Z.max(),
                            self.N.max(),
                            self.E.max(),
                        ]
                    ),
                )
                ylim = (
                    ylim1[0] + 0.1 * ylim1[0],
                    ylim1[1] + 0.1 * ylim1[1],
                )  # add 10% on each side to make it look nicer

        if jackshowall:
            subplots = True

        if subplots:
            if highf_tr is None:
                fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, figsize=(10, 8))
            else:
                fig, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=4, figsize=(10, 10))

                if infra_tr is not None:
                    plt.close(fig)
                    fig, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(
                        nrows=5, figsize=(10, 12)
                    )

            ax1.plot(tvec, self.Z, 'red', linewidth=1)
            ax1.set_ylabel('Up force (N)')
            ax2.plot(tvec, self.N, 'green', linewidth=1)
            ax2.set_ylabel('North force (N)')
            ax3.plot(tvec, self.E, 'blue', linewidth=1)
            ax3.set_ylabel('East force (N)')

            x = np.concatenate((tvec, tvec[::-1]))
            if self.jackknife is not None:
                if jackshowall:
                    for Z, N, E in zip(
                        self.jackknife.Z.all,
                        self.jackknife.N.all,
                        self.jackknife.E.all,
                    ):
                        ax1.plot(self.tvec, Z, 'red', alpha=0.2, linewidth=1)
                        ax2.plot(self.tvec, N, 'green', alpha=0.2, linewidth=1)
                        ax3.plot(self.tvec, E, 'blue', alpha=0.2, linewidth=1)
                else:
                    y = np.concatenate((Zlower, Zupper[::-1]))
                    poly = plt.Polygon(
                        list(zip(x, y)), facecolor='red', edgecolor='none', alpha=0.2
                    )
                    ax1.add_patch(poly)
                    y = np.concatenate((Nlower, Nupper[::-1]))
                    poly = plt.Polygon(
                        list(zip(x, y)), facecolor='green', edgecolor='none', alpha=0.2
                    )
                    ax2.add_patch(poly)
                    y = np.concatenate((Elower, Eupper[::-1]))
                    poly = plt.Polygon(
                        list(zip(x, y)), facecolor='blue', edgecolor='none', alpha=0.2
                    )
                    ax3.add_patch(poly)

            if highf_tr is not None:
                if type(highf_tr) != Trace:
                    raise TypeError('highf_tr is not an ObsPy Trace.')
                tvec2 = np.linspace(
                    0,
                    (len(highf_tr.data) - 1) * 1 / highf_tr.stats.sampling_rate,
                    num=len(highf_tr.data),
                )
                # Temporary fix, adjust for same zerotime
                if self.zero_time:
                    tvec2 -= self.zero_time
                tvec2 -= hfshift
                ax4.plot(tvec2, highf_tr.data * UMS_PER_MS, 'black')
                ax4.set_ylabel('Velocity (μm/s)')

            if infra_tr is not None:
                if type(infra_tr) != Trace:
                    raise TypeError('infra_tr is not an ObsPy Trace.')
                tvec2 = np.linspace(
                    0,
                    (len(infra_tr.data) - 1) * 1 / infra_tr.stats.sampling_rate,
                    num=len(infra_tr.data),
                )
                # Temporary fix, adjust for same zerotime
                if self.zero_time:
                    tvec2 -= self.zero_time
                tvec2 -= infra_shift
                ax5.plot(tvec2, infra_tr.data, 'black')
                ax5.set_ylabel('Pressure (Pa)')

                if infra_shift != 0:
                    ax5.annotate(
                        f'{infra_tr.id} (shifted –{infra_shift:1.0f} s)',
                        **annot_kwargs,
                    )
                else:
                    ax5.annotate(infra_tr.id, **annot_kwargs)

                # Remove x-axis labels for plots above this one
                for axis in [ax1, ax2, ax3, ax4]:
                    plt.setp(axis.get_xticklabels(), visible=False)

            if not xlim:
                xlim = [self.tvec.min(), self.tvec.max()]
            axes = fig.get_axes()
            [axe.set_xlim(xlim) for axe in axes]
            [axe.grid(True) for axe in axes]

            if sameY or ylim is not None:
                ax1.set_ylim(ylim)
                ax2.set_ylim(ylim)
                ax3.set_ylim(ylim)

            if self.impose_zero:
                [axe.axvline(0, color='gray', linestyle='solid', lw=3) for axe in axes]
            if self.max_duration is not None:
                [
                    axe.axvline(
                        self.max_duration, color='gray', linestyle='solid', lw=3
                    )
                    for axe in axes
                ]
            if highf_tr is not None:
                if hfylabel is not None:
                    ax4.set_ylabel(hfylabel)
                if hfshift != 0:
                    ax4.annotate(
                        f'{highf_tr.id} (shifted –{hfshift:1.0f} s)', **annot_kwargs
                    )
                else:
                    ax4.annotate(highf_tr.id, **annot_kwargs)

        else:
            if highf_tr is None:
                fig, ax = plt.subplots(figsize=(10, 4))
            else:
                fig, (ax, ax4) = plt.subplots(nrows=2, figsize=(10, 6))
                if type(highf_tr) != Trace:
                    raise TypeError('highf_tr is not an ObsPy Trace.')
                tvec2 = np.linspace(
                    0,
                    (len(highf_tr.data) - 1) * 1 / highf_tr.stats.sampling_rate,
                    num=len(highf_tr.data),
                )
                # Temporary fix, adjust for same zerotime
                if self.zero_time:
                    tvec2 -= self.zero_time
                tvec2 -= hfshift
                ax4.plot(tvec2, highf_tr.data * UMS_PER_MS, 'black')
                ax4.set_ylabel('Velocity (μm/s)')
                if hfshift != 0:
                    ax4.annotate(
                        f'{highf_tr.id} - shifted –{hfshift:1.0f} s', **annot_kwargs
                    )
                else:
                    ax4.annotate(highf_tr.id, **annot_kwargs)

            ax.plot(tvec, self.Z, 'red', label='Up')
            ax.plot(tvec, self.N, 'green', label='North')
            ax.plot(tvec, self.E, 'blue', label='East')

            if self.jackknife is not None:
                x = np.concatenate((tvec, tvec[::-1]))

                y = np.concatenate((Zlower, Zupper[::-1]))
                poly = plt.Polygon(
                    list(zip(x, y)), facecolor='red', edgecolor='none', alpha=0.2
                )
                ax.add_patch(poly)
                y = np.concatenate((Nlower, Nupper[::-1]))
                poly = plt.Polygon(
                    list(zip(x, y)), facecolor='green', edgecolor='none', alpha=0.2
                )
                ax.add_patch(poly)
                y = np.concatenate((Elower, Eupper[::-1]))
                poly = plt.Polygon(
                    list(zip(x, y)), facecolor='blue', edgecolor='none', alpha=0.2
                )
                ax.add_patch(poly)
            if xlim:
                ax.set_xlim(xlim)

            if highf_tr is not None:
                ax4.set_xlim(ax.get_xlim())
                ax4.grid(True)
            if hfylabel is not None:
                ax4.set_ylabel(hfylabel)

            ax.legend(loc='upper right')
            ax.grid(True)
            ax.set_ylabel('Force (N)')
            ax.set_ylim(ylim)
            if self.impose_zero:
                ax.axvline(0, color='gray', linestyle='solid', lw=3)
            if self.max_duration is not None:
                ax.axvline(self.max_duration, color='gray', linestyle='solid', lw=3)

        t0 = self.data.st_proc[0].stats.starttime
        if self.zero_time:
            t0 += self.zero_time
        plt.xlabel('Time (s) from {}'.format(t0.strftime('%Y-%m-%d %H:%M:%S')))

        fig.tight_layout()
        fig.show()

        return fig

    def plot_angle_magnitude(self, xlim=None, ylim=None, tvecshift=0.0):
        r"""Plot angles and magnitudes of inversion result.

        Args:
            xlim (list or tuple): x-axis limits
            ylim (list or tuple): y-axis limits
            tvecshift (int or float): [s] Shift time vector manually by this many
                seconds, not usually needed, only for display purposes (REMOVE?)

        Returns:
            :class:`~matplotlib.figure.Figure`: Output figure handle
        """

        tvec = self.tvec - tvecshift

        if self.jackknife is None:
            if ylim is None:
                ylim1 = (
                    np.amin([self.Z.min(), self.E.min(), self.N.min()]),
                    np.amax([self.Z.max(), self.E.max(), self.N.max()]),
                )
                ylim = (
                    ylim1[0] + 0.1 * ylim1[0],
                    ylim1[1] + 0.1 * ylim1[1],
                )  # add 10% on each side to make it look nicer
        else:
            Zupper = self.jackknife.Z.upper
            Nupper = self.jackknife.N.upper
            Eupper = self.jackknife.E.upper
            Zlower = self.jackknife.Z.lower
            Nlower = self.jackknife.N.lower
            Elower = self.jackknife.E.lower

            if ylim is None:
                ylim1 = (
                    np.amin([Zlower.min(), Elower.min(), Nlower.min()]),
                    np.amax([Zupper.max(), Eupper.max(), Nupper.max()]),
                )
            ylim = (
                ylim1[0] + 0.1 * ylim1[0],
                ylim1[1] + 0.1 * ylim1[1],
            )  # add 10% on each side to make it look nicer
        fig = plt.figure(figsize=(10, 10))

        # Plot the inversion result in the first one
        ax = fig.add_subplot(411)
        ax.plot(tvec, self.Z, 'red', label='Up')
        ax.plot(tvec, self.N, 'green', label='North')
        ax.plot(tvec, self.E, 'blue', label='East')

        if self.jackknife is not None:
            x = np.concatenate((tvec, tvec[::-1]))
            y = np.concatenate((Zlower, Zupper[::-1]))
            poly = plt.Polygon(
                list(zip(x, y)), facecolor='red', edgecolor='none', alpha=0.2
            )
            ax.add_patch(poly)
            y = np.concatenate((Nlower, Nupper[::-1]))
            poly = plt.Polygon(
                list(zip(x, y)), facecolor='green', edgecolor='none', alpha=0.2
            )
            ax.add_patch(poly)
            y = np.concatenate((Elower, Eupper[::-1]))
            poly = plt.Polygon(
                list(zip(x, y)), facecolor='blue', edgecolor='none', alpha=0.2
            )
            ax.add_patch(poly)

        if xlim:
            ax.set_xlim(xlim)

        ax.legend(loc='upper right')
        ax.grid(True)
        ax.set_ylabel('Force (N)')
        ax.set_ylim(ylim)

        # Plot the magnitudes in second one
        ax1 = fig.add_subplot(412)
        Mag = np.linalg.norm(list(zip(self.Z, self.E, self.N)), axis=1)
        ax1.plot(tvec, Mag, 'k', label='best')

        if self.jackknife is not None:
            MagU = np.linalg.norm(
                list(
                    zip(
                        np.maximum(np.abs(Zupper), np.abs(Zlower)),
                        np.maximum(np.abs(Eupper), np.abs(Elower)),
                        np.maximum(np.abs(Nupper), np.abs(Nlower)),
                    )
                ),
                axis=1,
            )
            MagL = np.linalg.norm(
                list(
                    zip(
                        np.minimum(np.abs(Zupper), np.abs(Zlower)),
                        np.minimum(np.abs(Eupper), np.abs(Elower)),
                        np.minimum(np.abs(Nupper), np.abs(Nlower)),
                    )
                ),
                axis=1,
            )
            ax1.plot(tvec, MagL, 'red', label='lower')
            ax1.plot(tvec, MagU, 'red', label='upper')
        else:
            MagU = None
            MagL = None
        ax1.set_ylabel('Force (N)')
        ax1.set_ylim(bottom=0)

        # Plot the horizontal azimuth
        ax2 = fig.add_subplot(413)
        tempang = (180 / np.pi) * np.arctan2(
            self.N, self.E
        ) - 90  # get angle counterclockwise relative to N
        # any negative values, add 360
        for i, temp in enumerate(tempang):
            if temp < 0:
                tempang[i] = temp + 360
        if self.jackknife is not None:
            tempangU = (180 / np.pi) * np.arctan2(Nupper, Eupper) - 90
            for i, temp in enumerate(tempangU):
                if temp < 0:
                    tempangU[i] = temp + 360
            tempangL = (180 / np.pi) * np.arctan2(Nlower, Elower) - 90
            for i, temp in enumerate(tempangL):
                if temp < 0:
                    tempangL[i] = temp + 360
        # now flip to clockwise to get azimuth
        Haz = 360 - tempang
        ax2.plot(tvec, Haz)
        ax2.set_ylabel('Azimuth (deg CW from N)')
        ax2.set_ylim(0, 360)

        # Plot the vertical angle
        ax3 = fig.add_subplot(414)
        Vang = (180 / np.pi) * np.arctan(self.Z / np.sqrt(self.N ** 2 + self.E ** 2))
        ax3.plot(tvec, Vang)
        ax3.set_ylabel('Vertical angle (deg)')

        axes = fig.get_axes()
        if xlim:
            axes = fig.get_axes()
            [axe.set_xlim(xlim) for axe in axes]
            [axe.grid(True) for axe in axes]

        if self.impose_zero:
            [axe.axvline(0, color='gray', linestyle='solid', lw=3) for axe in axes]
        if self.max_duration is not None:
            [
                axe.axvline(self.max_duration, color='gray', linestyle='solid', lw=3)
                for axe in axes
            ]

        plt.xlabel('Time (s)')

        self.angle_magnitude = AttribDict(
            magnitude=Mag,
            magnitude_upper=MagU,
            magnitude_lower=MagL,
            vertical_angle=Vang,
            horizontal_angle=Haz,
        )

        fig.tight_layout()
        fig.show()

        return fig

    def saverun(
        self,
        filepath=None,
        timestamp=False,
        figs2save=None,
        figs2save_names=None,
        light=True,
        filetype='png',
    ):
        r"""Save a force inversion run for later use.

        Args:
            filepath (str): Full path to directory where all files should be saved. If
                `None`, will use `self.main_folder`
            timestamp (bool): Name results with current time to avoid overwriting
                previous results
            figs2save (list or tuple): Figure handles to save
            figs2save_names (list or tuple): Names of figures (appends to end)
            light (bool): If `True`, does not save seismic data with object to save size
            filetype (str): Filetype given as extension, e.g. `'png'`
        """

        if filepath is None:
            filepath = self.main_folder

        obj = copy.copy(self)
        if light:
            obj.st = None

        # Get file name prefix
        if self.jackknife is None:
            jk = ''
        else:
            jk = 'JK'

        if timestamp:
            filename = '{}_{:1.0f}-{:1.0f}sec_{}{}'.format(
                self.nickname,
                self.filter['periodmin'],
                self.filter['periodmax'],
                jk,
                UTCDateTime.now().strftime('%Y-%m-%dT%H%M'),
            )
        else:
            filename = '{}_{:1.0f}-{:1.0f}sec_{}'.format(
                self.nickname, self.filter['periodmin'], self.filter['periodmax'], jk,
            )

        with open(os.path.join(filepath, f'{filename}.pickle'), 'wb') as f:
            pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

        if figs2save is not None:
            if figs2save_names is None:
                figs2save_names = range(len(figs2save))
            for i, fig in enumerate(figs2save):
                fig.savefig(
                    os.path.join(
                        filepath, f'{filename}_{figs2save_names[i]}.{filetype}'
                    ),
                    bbox_inches='tight',
                )


def _find_alpha(
    Ghat, dhat, I, L1=0, L2=0, tikhonov_ratios=(1.0, 0.0, 0.0), rough=False
):
    r"""Finds best regularization (trade-off) parameter alpha.

    Computes model with many values of alpha, plots L-curve, and finds point of steepest
    curvature where slope is negative.

    Args:
        Ghat (array): (m x n) matrix
        dhat (array): (1 x n) array of weighted data
        I (array): Identity matrix
        L1 (array): First order roughening matrix. If `0`, will use only 0th-order
            Tikhonov regularization
        L2 (array): Second order roughening matrix. If `0`, will use only 0th-order
            Tikhonov regularization
        tikhonov_ratios (list or tuple): Proportion each regularization method
            contributes to the overall regularization effect, where values correspond to
            [0th order, 1st order, 2nd order]. Must sum to 1
        rough (bool): If `False`, will do two iterations to fine tune the alpha
            parameter. If `True`, time will be saved because it will only do one round
            of searching

    Returns:
        tuple: Tuple containing:

        - **bestalpha** (float) – The optimal alpha
        - **fit1** (1D array) – List of residuals
        - **size1** (1D array) – List of model norms
        - **alphas** (1D array) – List of alphas tried
    """

    # Roughly estimate largest singular value (should not use alpha larger than expected
    # largest singular value)
    templ1 = np.ceil(np.log10(np.linalg.norm(Ghat)))
    templ2 = np.arange(templ1 - 5, templ1 - 1, 0.5)
    alphas = 10 ** templ2
    fit1 = []
    size1 = []

    Apart = Ghat.conj().T @ Ghat
    if type(L2) == float or type(L2) == int:
        L2part = 0.0
    else:
        L2part = L2.T @ L2
    if type(L1) == float or type(L1) == int:
        L1part = 0.0
    else:
        L1part = L1.T @ L1

    x = np.squeeze(Ghat.conj().T @ dhat)

    if rough:
        maxloops = 1
    else:
        maxloops = 2
    loop = 1
    while loop <= maxloops:
        for alpha in alphas:
            A = Apart + alpha ** 2 * (
                tikhonov_ratios[0] * I
                + tikhonov_ratios[1] * L1part
                + tikhonov_ratios[2] * L2part
            )  # Combo of all regularization things
            model, residuals, rank, s = sp.linalg.lstsq(A, x)

            temp1 = Ghat @ model.T - dhat
            fit1.append(sp.linalg.norm(temp1))
            size1.append(
                sp.linalg.norm(tikhonov_ratios[0] * model)
                + sp.linalg.norm(tikhonov_ratios[1] * np.dot(L1part, model))
                + sp.linalg.norm(tikhonov_ratios[2] * np.dot(L2part, model))
            )
        fit1 = np.array(fit1)
        size1 = np.array(size1)

        curves = _curvature(np.log10(fit1), np.log10(size1))
        # Zero out any points where function is concave to avoid picking points from
        # drop-off at end
        slp2 = np.gradient(np.gradient(np.log10(size1), np.log10(fit1)), np.log10(fit1))
        alphas = np.array(alphas)
        tempcurve = curves.copy()
        tempcurve[slp2 < 0] = np.inf
        idx = np.argmin(tempcurve)
        bestalpha = alphas[idx]
        loop += 1
        if loop > maxloops:
            break
        else:  # Loop again over smaller range
            alphas = np.logspace(
                np.round(np.log10(bestalpha)) - 1, np.round(np.log10(bestalpha)) + 1, 12
            )
            fit1 = []
            size1 = []

    _Lcurve(fit1, size1, alphas, bestalpha=bestalpha)
    if type(bestalpha) == list:
        if len(bestalpha) > 1:
            raise ValueError('Returned more than one alpha value, check codes.')
        bestalpha = bestalpha[0]

    return bestalpha, fit1, size1, alphas


def _Lcurve(fit1, size1, alphas, bestalpha=None):
    r"""Plot an L-curve.

    Args:
        fit1 (1D array): List of residuals
        size1 (1D array): List of model norms
        alphas (1D array): List of alphas tried
        bestalpha (float): The alpha value chosen
    """

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.loglog(fit1, size1, '.', markersize=10, color='black')
    for i, alpha in enumerate(alphas):
        ax.text(fit1[i], size1[i], f'  {alpha:.1e}', va='center')
    if bestalpha is not None:
        idx = np.argmin(np.abs(alphas - bestalpha))
        ax.plot(fit1[idx], size1[idx], 'or')
    ax.set_xlabel(r'Residual norm $||{\bfG}{\bfm}-{\bfd}||^2$')
    ax.set_ylabel(r'Solution norm $||{\bfm}||^2$')

    fig.tight_layout()
    fig.show()


def _varred(dt, dtnew):
    r"""Compute variance reduction :math:`\mathrm{VR}`.

    The formula is

    .. math::

        \mathrm{VR} = \left(1 - \frac{\|\mathbf{d}
        - \mathbf{d}_\mathbf{obs}\|^2}{\|\mathbf{d}_\mathbf{obs}\|^2}\right)
        \times 100\%\,,

    where :math:`\mathbf{d}_\mathbf{obs}` are the observed data, `dt`, and
    :math:`\mathbf{d}` are the synthetic data predicted by the forward model, `dtnew`.

    Args:
        dt: Array of original data
        dtnew: Array of modeled data

    Returns:
        float: Variance reduction :math:`\mathrm{VR}`
    """

    shp = np.shape(dt)
    shp = shp[0] * shp[1]
    dt_temp = np.reshape(dt, shp)
    dtnew_temp = np.reshape(dtnew, shp)
    d_dnew2 = (dt_temp - dtnew_temp) ** 2
    d2 = dt_temp ** 2
    VR = (1 - (np.sum(d_dnew2) / np.sum(d2))) * 100

    return VR


def _makeconvmat(c, size=None):
    r"""Build matrix used for convolution as implemented by matrix multiplication.

    Args:
        c (1D array): Signal to make convolution matrix for
        size (list or tuple): Optional input for desired size as ``(nrows, ncols)``;
            this will just shift ``cflip`` until it reaches the right size

    Returns:
        :class:`~numpy.ndarray`: Convolution matrix
    """

    cflip = c[::-1]  # flip order
    if size is None:
        C = np.zeros((2 * len(c) - 1, 2 * len(c) - 1))
        for i in range(2 * len(c) - 1):
            if i > len(c) - 1:
                zros = np.zeros(i + 1 - len(c))
                p = np.concatenate((zros, cflip, np.zeros(2 * len(c))))
            else:
                p = np.concatenate(((cflip[-(i + 1) :]), np.zeros(2 * len(c))))
            p = p[: 2 * len(c) - 1]
            C[i, :] = p.copy()
    else:
        # make it the right size
        C = np.zeros(size)
        for i in range(size[0]):
            if i > len(c) - 1:
                zros = np.zeros(i + 1 - len(c))
                p = np.concatenate((zros, cflip, np.zeros(size[1])))
            else:
                p = np.concatenate(((cflip[-(i + 1) :]), np.zeros(size[1])))
            p = p[: size[1]]  # cut p to the right size
            C[i, :] = p.copy()

    return C


def _makeshiftmat(c, shiftby, size1):
    r"""Build matrix that can be used for shifting of overlapping triangles.

    Used for triangle method. Signal goes across rows and each shift is a new column
    (opposite orientation to :func:`_makeconvmat()`)

    Args:
        c: Array of data (usually Green's function)
        shiftby (int): Number of samples to shift Green's function in each row
        size1 (list or tuple): Shape ``(nrows, ncols)`` of desired result. Will pad `c`
            if ``nrows`` is greater than ``len(c)``. Will shift `c` forward `shiftby`
            :math:`\times` ``ncols`` times

    Returns:
        :class:`~numpy.ndarray`: Matrix of shifted `c` of size `size1`
    """

    diff = len(c) - size1[0]
    if diff < 0:
        cpad = np.pad(c, (0, -diff), mode='edge')
    elif diff > 0:
        cpad = c[: size1[0]]
    else:
        cpad = c
    C = np.zeros(size1)
    for i in range(size1[1]):  # Loop over shifts and apply
        nshift = i * shiftby
        temp = np.pad(cpad.copy(), (nshift, 0), mode='edge')  # , end_values=(0., 0.))
        temp = temp[: size1[0]]
        C[:, i] = temp.copy()

    return C


def _curvature(x, y):
    r"""Estimate radius of curvature for each point on line to find corner of L-curve.

    Args:
        x: Array of x data
        y: Array of y data

    Returns:
        :class:`~numpy.ndarray`: Radius of curvature for each point (ends will be NaN)
    """

    # FOR EACH SET OF THREE POINTS, FIND RADIUS OF CIRCLE THAT FITS THEM - IGNORE ENDS
    R_2 = np.ones(len(x)) * float(
        'inf'
    )  # end numbers should be infinity because is straight line
    for i in range(1, len(R_2) - 1):
        xsub = x[i - 1 : i + 2]
        ysub = y[i - 1 : i + 2]
        m1 = -1 / (
            (ysub[0] - ysub[1]) / (xsub[0] - xsub[1])
        )  # slope of bisector of first segment
        m2 = -1 / (
            (ysub[1] - ysub[2]) / (xsub[1] - xsub[2])
        )  # slope of bisector of second segment
        b1 = ((ysub[0] + ysub[1]) / 2) - m1 * (
            (xsub[0] + xsub[1]) / 2
        )  # compute b for first bisector
        b2 = ((ysub[1] + ysub[2]) / 2) - m2 * (
            (xsub[1] + xsub[2]) / 2
        )  # compute b for second bisector

        Xc = (b1 - b2) / (m2 - m1)  # find intercept point of bisectors
        Yc = b2 + m2 * Xc

        R_2[i] = np.sqrt(
            (xsub[0] - Xc) ** 2 + (ysub[0] - Yc) ** 2
        )  # get distance from any point to intercept of bisectors to get radius

    return R_2


def _read(url):
    r"""Wrapper for :func:`obspy.core.stream.read` for long URLs."""
    with tempfile.NamedTemporaryFile() as f:
        urlretrieve(url, f.name)
        return read(f.name)
