import copy
import glob
import math
import os
import pickle
import random as rnd
import shutil
import stat
import subprocess

import numpy as np
import scipy as sp
from matplotlib import lines as mlines
from matplotlib import pyplot as plt
from obspy import Trace, UTCDateTime, read
from obspy.core import AttribDict
from obspy.signal.util import next_pow_2


class LSForce:
    r"""Class for performing force inversions.

    TODO:
        Make sure all attributes are defined in :meth:`__init__()` and complete this
        docstring!

    Attributes:
        data (:class:`~lsforce.lsdata.LSData`): The LSData object associated with this
            inversion
        domain: TODO
        data_sampling_rate: TODO
        nickname: TODO
        gf_computed: TODO
        gf_length: TODO
        inversion_complete: TODO
        main_folder: TODO
        method: TODO
        model_file: TODO
        T0: TODO
        triangle_half_width: TODO
        gf_sac_dir: TODO
        gf_run_dir: TODO
        filter: TODO
        data_length: TODO
        force_sampling_rate: TODO
        W: TODO
        Wvec: TODO
        weights: TODO
        add_to_zero: TODO
        zero_time: TODO
        impose_zero: TODO
        max_duration: TODO
        jackknife: TODO
        angle_magnitude: TODO
        G: TODO
        d: TODO
        model: Model vector of concatenated components (n x 1) of solution
        Z: [N] Vertical force time series extracted from model (positive up)
        N: [N] North force time series extracted from model (positive north)
        E: [N] East force time series extracted from model (positive east)
        tvec: [s] Time vector, referenced using `zero_time` (if specified) and corrected
            for `T0` time shift
        VR: [%] Variance reduction. Rule of thumb: This should be ~50–80%, if ~100%,
            solution is fitting data exactly and results are suspect. If ~5%, model may
            be wrong or something else may be wrong with setup
        dtorig: Original data vector (time domain)
        dtnew: Modeled data vector (Gm-d) (converted to time domain if domain is
            `'frequency'`)
        alpha: Regularization parameter that was used
        alphafit: TODO
    """

    def __init__(
        self,
        data,
        data_sampling_rate,
        domain='time',
        nickname=None,
        main_folder=None,
        method='tik',
    ):
        r"""Create an LSForce object.

        TODO:
            Implement sinusoid method!

        Args:
            data (:class:`~lsforce.lsdata.LSData`): LSData object, corrected for station
                response but not filtered
            data_sampling_rate (int or float): [Hz] Samples per second to use in
                inversion. All data will be resampled to this rate, and Green's
                functions will be created with this rate
            domain (str): Domain in which to do inversion, one of `'time'` or
                `'frequency'`
            nickname (str): Nickname for this event, used for convenient naming of files
            main_folder (str): If `None`, will use current folder
            method (str): One of `'tik'` — full waveform inversion using Tikhonov
                regularization (L2 norm minimization); `'triangle'` — inversion
                parameterized using overlapping triangles, variation on method of
                Ekström & Stark (2013); `'basis'` — parameterized using many Hanning
                basis functions; `'sinusoid'` — parameterized using a single sinusoid,
                variation on method of ????
        """

        self.data = data
        self.domain = domain
        self.data_sampling_rate = data_sampling_rate
        self.nickname = nickname
        self.gf_computed = False
        self.inversion_complete = False

        if main_folder is None:
            self.main_folder = os.getcwd()
        else:
            self.main_folder = main_folder

        if method not in ['tik', 'triangle']:
            raise ValueError(f'Method {method} not yet implemented.')

        self.method = method
        if self.method == 'triangle' and self.domain == 'frequency':
            raise ValueError('The triangle method must be done in the time domain.')

    def compute_greens(self, model_file, gf_duration, T0, triangle_half_width=5.0):
        r"""Compute Green's functions for inversion.

        Use CPS to compute the necessary Green's functions (GFs) for the source location
        and seismic stations being used. Computes the type of GFs appropriate for the
        method defined during class creation.

        Args:
            model_file (str): Full path to location of CPS model file
            gf_duration (int or float): [s] Duration of GFs
            T0 (int or float): [s] Amount of extra time prior to impulse application
            triangle_half_width (int or float): [s] Half width of isosceles triangle.
                Only needed for triangle method. This relates to the sampling interval
                of the force-time function since the triangles overlap by 50%
        """

        self.model_file = model_file
        self.T0 = T0

        self.triangle_half_width = triangle_half_width

        if self.nickname is None:
            self.nickname = ''

        self.gf_run_dir = os.path.join(
            self.main_folder,
            f'{self.nickname}_{os.path.splitext(os.path.basename(model_file))[0]}',
        )
        tmp_sac_dir = os.path.join(self.gf_run_dir, f'sacorig_{self.method}')
        self.gf_sac_dir = os.path.join(self.gf_run_dir, f'sacdata_{self.method}')

        # Make all the directories
        if not os.path.exists(self.main_folder):
            os.mkdir(self.main_folder)
        if not os.path.exists(self.gf_run_dir):
            os.mkdir(self.gf_run_dir)
        if not os.path.exists(tmp_sac_dir):
            os.mkdir(tmp_sac_dir)
        if not os.path.exists(self.gf_sac_dir):
            os.mkdir(self.gf_sac_dir)

        # write T0 file
        with open(os.path.join(self.gf_run_dir, 'T0.txt'), 'w') as f:
            f.write(f'{T0:3.2f}')

        # write triangle_half_width file, if applicable
        if self.method == 'triangle':
            with open(
                os.path.join(self.gf_run_dir, 'triangle_half_width.txt'), 'w'
            ) as f:
                f.write(f'{triangle_half_width:3.2f}')

        # Make sure there is only one occurrence of each station in list (ignore
        # channels)
        stacods = np.unique([tr.stats.station for tr in self.data.st_proc])
        dists = [
            self.data.st_proc.select(station=sta)[0].stats.distance for sta in stacods
        ]

        # write stadistlist.txt
        f = open(os.path.join(self.gf_run_dir, 'stadistlist.txt'), 'w')
        for sta, dis in zip(stacods, dists):
            f.write(f'{sta}\t{dis:5.1f}\n')
        f.close()

        # write dist file in free format
        # figure out how many samples
        print(f'Requested GF length = {gf_duration:g} s')
        samples = next_pow_2(gf_duration * self.data_sampling_rate)
        print(f'Optimized GF length = {samples / self.data_sampling_rate:g} s')
        f = open(os.path.join(self.gf_run_dir, 'dist'), 'w')
        for dis in dists:
            f.write(f'{dis:0.1f} {self.data_sampling_rate:0.2f} {samples:d} {T0:d} 0\n')
        f.close()
        self.gf_length = samples

        # move copy of model_file to current directory for recordkeeping
        shutil.copy2(
            model_file, os.path.join(self.gf_run_dir, os.path.basename(model_file))
        )

        # write shell script to run Green's functions
        shellscript = os.path.join(self.gf_run_dir, 'CPScommands.sh')
        with open(shellscript, 'w') as f:
            f.write('#!/bin/bash\n')
            f.write(f'rm {os.path.join(tmp_sac_dir, "*.sac")}\n')
            f.write(f'rm {os.path.join(self.gf_sac_dir, "*.sac")}\n')
            f.write(f'hprep96 -HR 0. -HS 0. -M {self.model_file} -d dist -R -EXF\n')
            f.write('hspec96 > hspec96.out\n')
            # 10^15 dynes is the default pulse here, equal to 10^10 Newtons
            if self.method == 'triangle':
                f.write(
                    'hpulse96 -d dist -D -t -l {} > Green\n'.format(
                        int(self.triangle_half_width / self.data_sampling_rate)
                    )
                )
            else:
                f.write('hpulse96 -d dist -V -OD -p > Green\n')
            f.write('f96tosac Green\n')
            f.write(f'cp *.sac {os.path.join(self.gf_sac_dir, ".")}\n')
            f.write(f'mv *.sac {os.path.join(tmp_sac_dir, ".")}\n')

        os.chmod(shellscript, stat.S_IRWXU)

        # Now actually run the codes
        currentdir = os.getcwd()
        os.chdir(self.gf_run_dir)
        proc = subprocess.Popen(
            shellscript, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        stdout, stderr = proc.communicate()
        retcode = proc.returncode
        if retcode != 0:
            os.chdir(currentdir)  # Change back to previous directory
            print(stderr)
            raise OSError(f'Green\'s functions were not computed:\n{stderr}')

        # copy and rename files
        files = [
            os.path.basename(x)
            for x in glob.glob(os.path.join(self.gf_sac_dir, '*.sac'))
        ]
        files.sort()
        for file1 in files:
            # get number of event
            indx = int(file1[1:4]) - 1
            GFt = file1[6:9]
            # rename
            newname = (
                f'GF_{self.nickname}_{stacods[indx]}_{dists[indx]:1.0f}km_{GFt}.sac'
            )
            os.rename(
                os.path.join(self.gf_sac_dir, file1),
                os.path.join(self.gf_sac_dir, newname),
            )

        os.chdir(currentdir)
        self.gf_computed = True

    def load_greens(self, model_file):
        r"""Load Green's functions for inversion.

        If Green's functions (GFs) were already computed for this exact data selection,
        this simply loads them for setup.

        TODO:
            Add error catching in case the stations in st don't line up with computed
            GFs in folder!

        Args:
            model_file (str): Full path to location of CPS model file. Used to locate
                appropriate directory containing GFs
        """

        if self.nickname is None:
            self.nickname = ''

        self.gf_run_dir = os.path.join(
            self.main_folder,
            f'{self.nickname}_{os.path.splitext(os.path.basename(model_file))[0]}',
        )
        if os.path.exists(os.path.join(self.gf_run_dir, f'sacorig_{self.method}')):
            self.gf_sac_dir = os.path.join(self.gf_run_dir, f'sacdata_{self.method}')
        else:
            self.gf_sac_dir = os.path.join(self.gf_run_dir, 'sacdata')

        # read T0 file
        with open(os.path.join(self.gf_run_dir, 'T0.txt'), 'r') as f:
            self.T0 = float(f.read())

        if self.method == 'triangle':
            with open(
                os.path.join(self.gf_run_dir, 'triangle_half_width.txt'), 'r'
            ) as f:
                self.triangle_half_width = float(f.read())

        # Read a file to get gf_length
        temp = read(glob.glob(os.path.join(self.gf_sac_dir, '*RVF*.sac'))[0])
        self.gf_length = len(temp[0])
        self.gf_computed = True

    def setup(
        self,
        period_range,
        weights=None,
        noise_window_dur=None,
        filter_order=2,
        zerophase=False,
    ):
        r"""Loads in GFs and creates all necessary matrices.

        Args:
            period_range (list or tuple): [s] Bandpass filter corners
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

        # Create filter dictionary to keep track of filter used without
        # creating too many new attributes
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

        # Since GFs are computed for a 10^15 dyne impulse, this converts GFs to what
        # they'd be for a 1 dyne impulse
        K = 1e-15
        self.data_length = len(st[0].data)

        if self.gf_length > self.data_length:
            raise ValueError(
                'gf_length is greater than data_length. Reselect data and/or '
                'recompute Green\'s functions so that data is longer than Green\'s '
                'functions'
            )

        if self.domain == 'time':
            # TODO: ADD WAY TO ACCOUNT FOR WHEN GF_LENGTH IS LONGER THAN DATA_LENGTH -
            #  ACTUALLY SHOULD BE AS LONG AS BOTH ADDED TOGETHER TO AVOID WRAPPING ERROR
            lenUall = self.data_length * len(st)
        elif self.domain == 'frequency':
            # Needs to be the length of the two added together because convolution
            # length M+N-1
            nfft = next_pow_2(self.data_length)  # + gf_length)
            lenUall = nfft * len(st)
        else:
            raise ValueError(
                'domain not recognized. Must be \'time\' or \'frequency\'.'
            )

        if self.method == 'tik':

            self.force_sampling_rate = self.data_sampling_rate

            # initialize weighting matrices
            Wvec = np.ones(lenUall)
            indx = 0
            weight = np.ones(self.data.st_proc.count())

            n = self.data_length

            for i, trace in enumerate(st):
                # find component of st
                component = trace.stats.channel[2]
                station = trace.stats.station
                if component == 'Z':
                    zvf = read(os.path.join(self.gf_sac_dir, f'*_{station}_*ZVF.sac'))
                    if len(zvf) > 1:
                        raise ValueError(f'Found more than one ZVF GF for {station}.')
                    else:
                        zvf = zvf[0]
                    zhf = read(os.path.join(self.gf_sac_dir, f'*_{station}_*ZHF.sac'))
                    if len(zhf) > 1:
                        raise ValueError(f'Found more than one ZHF GF for {station}.')
                    else:
                        zhf = zhf[0]
                    # process the same way as st (except shouldn't need to resample)
                    zvf.detrend()
                    zvf.taper(max_percentage=0.05)
                    zvf.filter(
                        'bandpass',
                        freqmin=self.filter['freqmin'],
                        freqmax=self.filter['freqmax'],
                        corners=self.filter['order'],
                        zerophase=self.filter['zerophase'],
                    )

                    zhf.detrend()
                    zhf.taper(max_percentage=0.05)
                    zhf.filter(
                        'bandpass',
                        freqmin=self.filter['freqmin'],
                        freqmax=self.filter['freqmax'],
                        corners=self.filter['order'],
                        zerophase=self.filter['zerophase'],
                    )

                    if self.domain == 'time':
                        ZVF = _makeconvmat(zvf.data, size=(n, n))
                        ZHF = _makeconvmat(zhf.data, size=(n, n))
                    else:
                        zvff = np.fft.fft(zvf.data, nfft)
                        zhff = np.fft.fft(zhf.data, nfft)
                        ZVF = np.diag(zvff)
                        ZHF = np.diag(zhff)
                    az = math.radians(trace.stats.azimuth)
                    newline = np.hstack(
                        (K * ZVF, K * ZHF * math.cos(az), K * ZHF * math.sin(az))
                    )
                elif component == 'R':
                    rvf = read(os.path.join(self.gf_sac_dir, f'*_{station}_*RVF.sac'))
                    if len(rvf) > 1:
                        raise ValueError(f'Found more than one RVF GF for {station}.')
                    else:
                        rvf = rvf[0]
                    rhf = read(os.path.join(self.gf_sac_dir, f'*_{station}_*RHF.sac'))
                    if len(rhf) > 1:
                        raise ValueError(f'Found more than one RHF GF for {station}.')
                    else:
                        rhf = rhf[0]

                    # process the same way as st
                    rvf.detrend()
                    rvf.taper(max_percentage=0.05)

                    rvf.filter(
                        'bandpass',
                        freqmin=self.filter['freqmin'],
                        freqmax=self.filter['freqmax'],
                        corners=self.filter['order'],
                        zerophase=self.filter['zerophase'],
                    )

                    rhf.detrend()
                    rhf.taper(max_percentage=0.05)
                    rhf.filter(
                        'bandpass',
                        freqmin=self.filter['freqmin'],
                        freqmax=self.filter['freqmax'],
                        corners=self.filter['order'],
                        zerophase=self.filter['zerophase'],
                    )
                    if self.domain == 'time':
                        RVF = _makeconvmat(rvf.data, size=(n, n))
                        RHF = _makeconvmat(rhf.data, size=(n, n))
                    else:
                        rvff = np.fft.fft(rvf.data, nfft)
                        rhff = np.fft.fft(rhf.data, nfft)
                        RVF = np.diag(rvff)
                        RHF = np.diag(rhff)
                    az = math.radians(trace.stats.azimuth)
                    newline = np.hstack(
                        (K * RVF, K * RHF * math.cos(az), K * RHF * math.sin(az))
                    )
                elif component == 'T':
                    thf = read(os.path.join(self.gf_sac_dir, f'*_{station}_*THF.sac'))
                    if len(thf) > 1:
                        raise ValueError(f'Found more than one THF GF for {station}.')
                    else:
                        thf = thf[0]
                    # process the same way as st
                    thf.detrend()
                    thf.taper(max_percentage=0.05)
                    thf.filter(
                        'bandpass',
                        freqmin=self.filter['freqmin'],
                        freqmax=self.filter['freqmax'],
                        corners=self.filter['order'],
                        zerophase=self.filter['zerophase'],
                    )
                    if self.domain == 'time':
                        THF = _makeconvmat(thf.data, size=(n, n))
                    else:
                        thff = np.fft.fft(thf.data, nfft)
                        THF = np.diag(thff)
                    TVF = 0.0 * THF.copy()
                    az = math.radians(trace.stats.azimuth)
                    newline = np.hstack(
                        (TVF, K * THF * math.sin(az), -K * THF * math.cos(az))
                    )
                else:
                    raise ValueError(f'st not rotated to T and R for {station}.')
                # Deal with data
                if self.domain == 'time':
                    datline = trace.data
                else:
                    datline = np.fft.fft(trace.data, nfft)

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
                            trace.data[
                                0 : int(noise_window_dur * trace.stats.sampling_rate)
                            ]
                        )
                    elif weights == 'distance':
                        weight[i] = trace.stats.distance

                    Wvec[indx : indx + self.data_length] = (
                        Wvec[indx : indx + self.data_length] * weight[i]
                    )
                    indx += self.data_length

        elif self.method == 'triangle':

            # initialize weighting matrices
            Wvec = np.ones(lenUall)
            indx = 0
            weight = np.ones(self.data.st_proc.count())

            n = self.data_length
            fshiftby = int(
                self.triangle_half_width / self.data_sampling_rate
            )  # Number of samples to shift each triangle by
            Flen = int(
                np.floor(self.data_length / fshiftby)
            )  # Number of shifts, corresponds to length of force time function
            self.force_sampling_rate = 1.0 / fshiftby

            for i, trace in enumerate(st):
                # find component of st
                component = trace.stats.channel[2]
                station = trace.stats.station
                if component == 'Z':
                    zvf = read(os.path.join(self.gf_sac_dir, f'*{station}*ZVF.sac'))
                    if len(zvf) > 1:
                        raise ValueError(f'Found more than one ZVF GF for {station}.')
                    else:
                        zvf = zvf[0]
                    zhf = read(os.path.join(self.gf_sac_dir, f'*{station}*ZHF.sac'))
                    if len(zhf) > 1:
                        raise ValueError(f'Found more than one ZHF GF for {station}.')
                    else:
                        zhf = zhf[0]
                    # process the same way as st (except shouldn't need to resample)
                    """ Don't need to filter these GFs? Has non-zero offset so filtering
                    does weird things, already convolved with LP source-time function
                    zvf.detrend()
                    zvf.taper(max_percentage=0.05)
                    zvf.filter('bandpass', freqmin=self.filter['freqmin'],
                               freqmax=self.filter['freqmax'],
                               corners=self.filter['order'],
                               zerophase=self.filter['zerophase'])

                    zhf.detrend()
                    zhf.taper(max_percentage=0.05)
                    zhf.filter('bandpass', freqmin=self.filter['freqmin'],
                               freqmax=self.filter['freqmax'],
                               corners=self.filter['order'],
                               zerophase=self.filter['zerophase'])
                    """
                    ZVF = _makeshiftmat(zvf.data, shiftby=fshiftby, size1=(n, Flen))
                    ZHF = _makeshiftmat(zhf.data, shiftby=fshiftby, size1=(n, Flen))
                    az = math.radians(trace.stats.azimuth)
                    newline = np.hstack(
                        (K * ZVF, K * ZHF * math.cos(az), K * ZHF * math.sin(az))
                    )
                elif component == 'R':
                    rvf = read(os.path.join(self.gf_sac_dir, f'*{station}*RVF.sac'))
                    if len(rvf) > 1:
                        raise ValueError(f'Found more than one RVF GF for {station}.')
                    else:
                        rvf = rvf[0]
                    rhf = read(os.path.join(self.gf_sac_dir, f'*{station}*RHF.sac'))
                    if len(rhf) > 1:
                        raise ValueError(f'Found more than one RHF GF for {station}.')
                    else:
                        rhf = rhf[0]
                    """ Don't need to filter these GFs?
                    #process the same way as st
                    rvf.detrend()
                    rvf.taper(max_percentage=0.05)

                    rvf.filter('bandpass', freqmin=self.filter['freqmin'],
                               freqmax=self.filter['freqmax'],
                               corners=self.filter['order'],
                               zerophase=self.filter['zerophase'])

                    rhf.detrend()
                    rhf.taper(max_percentage=0.05)
                    rhf.filter('bandpass', freqmin=self.filter['freqmin'],
                               freqmax=self.filter['freqmax'],
                               corners=self.filter['order'],
                               zerophase=self.filter['zerophase'])
                    """
                    RVF = _makeshiftmat(rvf.data, shiftby=fshiftby, size1=(n, Flen))
                    RHF = _makeshiftmat(rhf.data, shiftby=fshiftby, size1=(n, Flen))
                    az = math.radians(trace.stats.azimuth)
                    newline = np.hstack(
                        (K * RVF, K * RHF * math.cos(az), K * RHF * math.sin(az))
                    )
                elif component == 'T':
                    thf = read(os.path.join(self.gf_sac_dir, f'*{station}*THF.sac'))
                    if len(thf) > 1:
                        raise ValueError(f'Found more than one THF GF for {station}.')
                    else:
                        thf = thf[0]
                    """ Don't need to filter these GFs?
                    #process the same way as st
                    thf.detrend()
                    thf.taper(max_percentage=0.05)
                    thf.filter('bandpass', freqmin=self.filter['freqmin'],
                               freqmax=self.filter['freqmax'],
                               corners=self.filter['order'],
                               zerophase=self.filter['zerophase'])
                    """
                    THF = _makeshiftmat(thf.data, shiftby=fshiftby, size1=(n, Flen))
                    TVF = 0.0 * THF.copy()
                    az = math.radians(trace.stats.azimuth)
                    newline = np.hstack(
                        (TVF, K * THF * math.sin(az), -K * THF * math.cos(az))
                    )
                else:
                    raise ValueError(f'st not rotated to T and R for {station}.')
                # Deal with data
                datline = trace.data

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
                            trace.data[
                                0 : int(noise_window_dur * trace.stats.sampling_rate)
                            ]
                        )
                    elif weights == 'distance':
                        weight[i] = trace.stats.distance

                    Wvec[indx : indx + self.data_length] = (
                        Wvec[indx : indx + self.data_length] * weight[i]
                    )
                    indx += self.data_length

        else:
            raise ValueError(f'Method {self.method} not supported.')

        # Normalize Wvec so largest weight is 1.
        self.Wvec = Wvec / np.max(np.abs(Wvec))
        self.weights = weight / np.max(np.abs(weight))

        if np.shape(G)[0] != len(d):
            raise ValueError(
                'G and d sizes are not compatible, fix something somewhere.'
            )
        # Need to multiply G by sample interval [s] since convolution is an integral
        self.G = G * 1.0 / self.data_sampling_rate
        self.d = d * 100.0  # Convert data from m to cm
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
        zero_scaler=15.0,
        zero_taper_length=20.0,
        tikhonov_ratios=(1.0, 0.0, 0.0),
    ):
        r"""Performs full-waveform inversion using Tikhonov regularization.

        Args:
            alphaset (int or float): Set regularization parameter. If `None`, will
                search for best alpha using the L-curve method
            zero_scaler (int or float): Factor by which to divide Gnorm to get scaling
                factor used for zero constraint. The lower the number, the stronger the
                constraint, but the higher the risk of high frequency oscillations due
                to a sudden release of the constraint
            zero_taper_length (int or float): [s] Length of taper for `zero_scaler`.
                Shorter tapers can result in sharp artifacts, so longer is better
            tikhonov_ratios (list or tuple): Proportion each regularization method
                contributes to the overall regularization effect, where values
                correspond to [0th order, 1st order, 2nd order]. Must sum to 1
        """

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
        gl = int(n / 3)  # self.data_length

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

        scaler = Ghatnorm / zero_scaler
        if self.impose_zero:  # tell model when there should be no forces
            # TODO get this to work for triangle method (need to change len methods)
            len2 = int(
                np.floor(((self.zero_time + self.T0) * self.force_sampling_rate))
            )
            if self.method == 'triangle':
                len2 = int(
                    np.floor(((self.zero_time + self.T0) * self.force_sampling_rate))
                )
            if self.method == 'tik':
                len3 = int(
                    zero_taper_length * self.force_sampling_rate
                )  # make it constant
                temp = np.hanning(2 * len3)
                temp = temp[len3:]
                vals2 = np.hstack((np.ones(len2 - len3), temp))
            else:  # No taper
                vals2 = np.hstack((np.ones(len2), np.zeros(gl - len2)))

            for i, val in enumerate(vals2):
                first1 = np.zeros(3 * gl)
                second1 = first1.copy()
                third1 = first1.copy()
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
                Ghat, dhat, I, L1, L2, tikhonov_ratios=tikhonov_ratios, invmethod='lsq',
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

        if self.domain == 'frequency':
            model, residuals, rank, s = sp.linalg.lstsq(A, x)
            self.model = model.copy()
            div = len(model) / 3

            # Convert from dynes to newtons, flip so up is positive
            self.Z = -np.real(np.fft.ifft(model[0:div]) / 10 ** 5)
            self.N = np.real(np.fft.ifft(model[div : 2 * div]) / 10 ** 5)
            self.E = np.real(np.fft.ifft(model[2 * div :]) / 10 ** 5)

            # run forward model
            df_new = self.G @ model.T
            # convert d and df_new back to time domain
            dt, dtnew = _back2time(self.d, df_new, self.data.st_proc.count(), dl)
            self.dtorig = dt
            self.dtnew = dtnew

        else:  # domain is time
            model, residuals, rank, s = sp.linalg.lstsq(A, x)
            self.model = model.copy()
            div = int(len(model) / 3)

            # Convert from dynes to newtons, flip so up is positive
            self.Z = -model[0:div] / 10 ** 5
            self.N = model[div : 2 * div] / 10 ** 5
            self.E = model[2 * div :] / 10 ** 5

            dtnew = self.G.dot(model)
            self.dtnew = np.reshape(dtnew, (self.data.st_proc.count(), dl))
            self.dtorig = np.reshape(self.d, (self.data.st_proc.count(), dl))

        # compute variance reduction
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
        if self.zero_time is not None:
            tvec -= self.zero_time
        if self.method == 'triangle':
            # Shift so that peak of triangle function lines up with time of force
            # interval
            tvec += self.triangle_half_width
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

                if self.domain == 'frequency':
                    model, residuals, rank, s = sp.linalg.lstsq(Aj, xj)
                    div = len(model) / 3
                    Zf = -np.real(
                        np.fft.ifft(model[0:div]) / 10 ** 5
                    )  # convert from dynes to newtons, flip so up is positive
                    Nf = np.real(np.fft.ifft(model[div : 2 * div]) / 10 ** 5)
                    Ef = np.real(np.fft.ifft(model[2 * div :]) / 10 ** 5)
                    # run forward model
                    df_new = Gtemp @ model.T
                    # convert d and df_new back to time domain
                    dt, dtnew = _back2time(dtemp, df_new, numkeep, dl)

                else:  # domain is time
                    model, residuals, rank, s = sp.linalg.lstsq(Aj, xj)
                    div = int(len(model) / 3)
                    Zf = (
                        -model[0:div] / 10 ** 5
                    )  # convert from dynes to netwons, flip so up is positive
                    Nf = model[div : 2 * div] / 10 ** 5
                    Ef = model[2 * div :] / 10 ** 5
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
                ms2ums = 1e6
                ax4.plot(tvec2, highf_tr.data * ms2ums, 'black')
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
                ax4.plot(tvec2, highf_tr.data)
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
    Ghat,
    dhat,
    I,
    L1=0,
    L2=0,
    invmethod='lsq',
    tikhonov_ratios=(1.0, 0.0, 0.0),
    rough=False,
):
    r"""Finds best regularization (trade-off) parameter alpha.

    Computes model with many values of alpha, plots L-curve, and finds point of steepest
    curvature where slope is negative.

    TODO:
        Finish this docstring!

    Args:
        Ghat (array): (m x n) matrix
        dhat (array): (1 x n) array of weighted data
        I (array): Identity matrix
        L1 (array): First order roughening matrix. If `0`, will use only 0th-order
            Tikhonov regularization
        L2 (array): Second order roughening matrix. If `0`, will use only 0th-order
            Tikhonov regularization
        invmethod (str): `'lsq'` — use least squares (regular Tikhonov); `'nnls'` — use
            non-negative least squares
        tikhonov_ratios (list or tuple): Proportion each regularization method
            contributes to the overall regularization effect, where values correspond to
            [0th order, 1st order, 2nd order]. Must sum to 1
        rough (bool): If `False`, will do two iterations to fine tune the alpha
            parameter. If `True`, time will be saved because it will only do one round
            of searching

    Returns:
        tuple: Tuple containing:

        - **bestalpha** – TODO
        - **fit1** – TODO
        - **size1** – TODO
        - **alphas** – TODO
    """

    # Roughly estimate largest singular value (should not use alpha larger than expected
    # largest singular value)
    templ1 = np.ceil(np.log10(np.linalg.norm(Ghat)))
    templ2 = np.arange(templ1 - 6, templ1 - 2)
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

    # rough first iteration
    for alpha in alphas:
        A = Apart + alpha ** 2 * (
            tikhonov_ratios[0] * I
            + tikhonov_ratios[1] * L1part
            + tikhonov_ratios[2] * L2part
        )  # Combo of all regularization things
        if invmethod == 'lsq':
            model, residuals, rank, s = sp.linalg.lstsq(A, x)
        elif invmethod == 'nnls':
            model, residuals = sp.optimize.nnls(A, x)
        else:
            raise ValueError(f'Inversion method {invmethod} not recognized.')
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
    # Zero out any points where function is concave to avoid picking points from dropoff
    # at end
    slp2 = np.gradient(np.gradient(np.log10(size1), np.log10(fit1)), np.log10(fit1))
    alphas = np.array(alphas)
    tempcurve = curves.copy()
    tempcurve[slp2 < 0] = np.max(curves)
    idx = np.argmin(tempcurve)
    alpha = alphas[idx]

    if not rough:
        # Then hone in
        alphas = np.logspace(
            np.round(np.log10(alpha)) - 1, np.round(np.log10(alpha)) + 1, 10
        )
        fit1 = []
        size1 = []
        for newalpha in alphas:
            A = Apart + newalpha ** 2 * (
                tikhonov_ratios[0] * I
                + tikhonov_ratios[1] * L1part
                + tikhonov_ratios[2] * L2part
            )  # Combo of all regularization things
            if invmethod == 'lsq':
                model, residuals, rank, s = sp.linalg.lstsq(A, x)
            elif invmethod == 'nnls':
                model, residuals = sp.optimize.nnls(A, x)

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
        # dropoff at end
        slp2 = np.gradient(np.gradient(np.log10(size1), np.log10(fit1)), np.log10(fit1))
        alphas = np.array(alphas)
        tempcurve = curves.copy()
        tempcurve[slp2 < 0] = np.max(curves)
        idx = np.argmin(tempcurve)
        bestalpha = alphas[idx]
    else:
        bestalpha = alpha

    _Lcurve(fit1, size1, alphas)
    if type(bestalpha) == list:
        if len(bestalpha) > 1:
            raise ValueError('Returned more than one alpha value, check codes.')
        bestalpha = bestalpha[0]

    return bestalpha, fit1, size1, alphas


def _Lcurve(fit1, size1, alphas):
    r"""Plot an L-curve.

    TODO:
        Finish this docstring!

    Args:
        fit1: TODO
        size1: TODO
        alphas: TODO
    """

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.loglog(fit1, size1, '.', markersize=10, color='black')
    for i, alpha in enumerate(alphas):
        ax.text(fit1[i], size1[i], f'  {alpha:.1e}', va='center')
    ax.set_xlabel(r'Residual norm $||{\bfG}{\bfm}-{\bfd}||^2$')
    ax.set_ylabel(r'Solution norm $||{\bfm}||^2$')

    fig.tight_layout()
    fig.show()


def _varred(dt, dtnew):
    r"""Compute variance reduction, :math:`\mathrm{VR}`, in the time domain.

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


def _back2time(d, df_new, numsta, datlenorig):
    r"""Convert data back to the time domain and cut off zero padding.

    TODO:
        Finish this docstring!

    Args:
        d: TODO
        df_new: TODO
        numsta: TODO
        datlenorig: TODO

    Returns:
        tuple: Tuple containing:

        - **dt** – TODO
        - **dtnew** – TODO
    """

    datlength = int(len(d) / numsta)
    dfrsp = np.reshape(d, (numsta, datlength))
    dfnrsp = np.reshape(df_new, (numsta, datlength))
    dt = np.real(np.fft.ifft(dfrsp, axis=1))
    dt = dt[0:, 0:datlenorig]
    dtnew = np.real(np.fft.ifft(dfnrsp, axis=1))
    dtnew = dtnew[0:, 0:datlenorig]

    return dt, dtnew


def _makeconvmat(c, size=None):
    r"""Build matrix that used for convolution as implemented by matrix multiplication.

    `size` is optional input for desired size as ``(nrows, ncols)``; this will just
    shift ``cflip`` until it reaches the right size.

    TODO:
        Finish this docstring!

    Args:
        c: TODO
        size (list or tuple): TODO

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
