import numpy as np
from obspy import Trace, read, UTCDateTime
from obspy.signal.util import next_pow_2
import math
import glob
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import scipy as sp
import random as rnd
import pickle
import os
import stat
import shutil
import subprocess
import copy


class LSForce:
    """
    Class for single force inversions
    """

    def __init__(
        self,
        data,
        sampling_rate,
        domain='time',
        nickname=None,
        main_folder=None,
        method='tik',
    ):
        """
        Args:
            data (LSData): LSData object, corrected for station response but not filtered
            sampling_rate (float): Number of samples per second (Hz) to use in inversion.
                All data will be resampled to this rate and Greens functions
                will be created with this sample rate.
            domain (str): domain in which to do inversion, 'time' (default) or
                'freq'
            nickname (str): Nickname for this event, used for convenient in
                naming files
            main_folder (str): if None, will use current folder
            method (str): 'tik' = full waveform inversion using Tikhonov
                                regularization (L2 norm minimization)
                          'triangle' = parameterized inversion using overlapping
                                  triangles (variation of method of Ekstrom et al., 2013)
                          'basis' = parameterized using many hanning basis functions
                          'sinusoid' = parameterized using single sinusoid
                                  (variation of method by Chao et al. YEAR)
        """

        # General
        self.st = data.st_proc.copy()
        self.domain = domain
        self.sampling_rate = sampling_rate
        self.nickname = nickname
        self.numsta = len(self.st)
        self.greens_computed = False
        self.lat = data.source_lat
        self.lon = data.source_lon
        self.inversion_complete = False

        if main_folder is None:
            self.main_folder = os.getcwd()
        else:
            self.main_folder = main_folder

        if method not in ['tik', 'triangle']:
            raise ValueError(f'Method {method} not yet implemented.')

        self.method = method
        if self.method in ['triangle'] and self.domain == 'freq':
            raise ValueError('The triangle method must be done in the time domain.')

    def compute_greens(self, model_file, gf_duration, T0, L=5.0):
        """
        Use CPS to compute the necessary Greens functions relating the source
        location and the seismic stations in st
        Will compute the type of Greens functions appropriate for the method
        defined during class initiation

        Args:
            model_file (str): Full file path to location of CPS model file
            gf_duration (float):
            T0 (float): number of seconds prior to impulse application
            L (float): half width, in seconds, of triangle. Only needed for
                triangle method, default 5 sec. This will also correspond
                to sample interval of resulting force time function because
                triangles overlap by 50%
        """

        self.model_file = model_file
        self.T0 = T0

        self.L = L

        if self.nickname is None:
            self.nickname = ''

        self.moddir = os.path.join(
            self.main_folder,
            '%s_%s'
            % (self.nickname, os.path.splitext(os.path.basename(model_file))[0]),
        )
        self.sacodir = os.path.join(self.moddir, 'sacorig_%s' % self.method)
        self.sacdir = os.path.join(self.moddir, 'sacdata_%s' % self.method)

        # Make all the directories
        if not os.path.exists(self.main_folder):
            os.mkdir(self.main_folder)
        if not os.path.exists(self.moddir):
            os.mkdir(self.moddir)
        if not os.path.exists(self.sacodir):
            os.mkdir(self.sacodir)
        if not os.path.exists(self.sacdir):
            os.mkdir(self.sacdir)

        # write T0 file
        with open(os.path.join(self.moddir, 'T0.txt'), 'w') as f:
            f.write('%3.2f' % T0)

        # write L file, if applicable
        if self.method in ['triangle']:
            with open(os.path.join(self.moddir, 'L.txt'), 'w') as f:
                f.write('%3.2f' % L)

        # Make sure there is only one occurrence of each station in list (ignore channels)
        stacods = np.unique([tr.stats.station for tr in self.st])
        dists = [self.st.select(station=sta)[0].stats.distance for sta in stacods]

        # write stadistlist.txt
        f = open(os.path.join(self.moddir, 'stadistlist.txt'), 'w')
        for sta, dis in zip(stacods, dists):
            f.write('%s\t%5.1f\n' % (sta, dis))
        f.close()

        # write dist file in free format
        # figure out how many samples
        samples = next_pow_2(gf_duration * self.sampling_rate)
        f = open(os.path.join(self.moddir, 'dist'), 'w')
        for dis in dists:
            f.write(
                '%0.1f %0.2f %i %i 0\n'
                % (dis, 1.0 / self.sampling_rate, samples, self.T0)
            )
        f.close()
        self.greenlength = samples

        # move copy of model_file to current directory for recordkeeping
        shutil.copy2(
            model_file, os.path.join(self.moddir, os.path.basename(model_file))
        )

        # write shell script to run Green's functions
        self.shellscript = os.path.join(self.moddir, 'CPScommands.sh')
        with open(self.shellscript, 'w') as f:
            f.write('#!/bin/bash\n')
            f.write('rm %s\n' % os.path.join(self.sacodir, '*.sac'))
            f.write('rm %s\n' % os.path.join(self.sacdir, '*.sac'))
            f.write(
                'hprep96 -HR 0. -HS 0. -M %s -d %s -R -EXF\n'
                % (self.model_file, 'dist')
            )
            f.write('hspec96 > hspec96.out\n')
            if self.method == 'triangle':
                f.write(
                    'hpulse96 -d %s -V -D -t -l %d > Green\n'
                    % ('dist', int(self.L / self.sampling_rate))
                )
            else:
                f.write('hpulse96 -d %s -V -OD -p > Green\n' % 'dist')
            f.write('f96tosac Green\n')
            f.write('cp *.sac %s\n' % os.path.join(self.sacdir, '.'))
            f.write('mv *.sac %s\n' % os.path.join(self.sacodir, '.'))

        os.chmod(self.shellscript, stat.S_IRWXU)

        # Now actually run the codes
        currentdir = os.getcwd()
        os.chdir(self.moddir)
        proc = subprocess.Popen(
            self.shellscript, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        stdout, stderr = proc.communicate()
        retcode = proc.returncode
        if retcode != 0:
            os.chdir(currentdir)  # Change back to previous directory
            print(stderr)
            raise OSError(f'Green\'s functions were not computed:\n{stderr}')

        # copy and rename files
        files = [
            os.path.basename(x) for x in glob.glob(os.path.join(self.sacdir, '*.sac'))
        ]
        files.sort()
        for file1 in files:
            # get number of event
            indx = int(file1[1:4]) - 1
            GFt = file1[6:9]
            # rename
            newname = 'GF_%s_%s_%1.0fkm_%s.sac' % (
                self.nickname,
                stacods[indx],
                dists[indx],
                GFt,
            )
            os.rename(
                os.path.join(self.sacdir, file1), os.path.join(self.sacdir, newname)
            )

        os.chdir(currentdir)
        self.greens_computed = True

    def load_greens(self, model_file):

        """
        If Greens functions were already computed for this exact data
        selection, this simply loads info about them needed for setup

        Args:
            model_file (str): the name of the model file used to compute the
                Greens functions. This is so they can be found because they
                are saved in a folder referencing the model file name

            TODO add error catching in case the stations in st don't line
            up with computed GFs in folder
        """

        if self.nickname is None:
            self.nickname = ''

        self.moddir = os.path.join(
            self.main_folder,
            '%s_%s'
            % (self.nickname, os.path.splitext(os.path.basename(model_file))[0]),
        )
        if os.path.exists(os.path.join(self.moddir, 'sacorig_%s' % self.method)):
            self.sacodir = os.path.join(self.moddir, 'sacorig_%s' % self.method)
            self.sacdir = os.path.join(self.moddir, 'sacdata_%s' % self.method)
        else:
            self.sacodir = os.path.join(self.moddir, 'sacorig')
            self.sacdir = os.path.join(self.moddir, 'sacdata')

        # read T0 file
        with open(os.path.join(self.moddir, 'T0.txt'), 'r') as f:
            self.T0 = float(f.read())

        if self.method in ['triangle']:
            with open(os.path.join(self.moddir, 'L.txt'), 'r') as f:
                self.L = float(f.read())

        # Read a file to get greenlength
        temp = read(glob.glob(os.path.join(self.sacdir, '*RVF*.sac'))[0])
        self.greenlength = len(temp[0])
        self.greens_computed = True

    def setup(
        self,
        weights=None,
        weightpre=None,
        period_range=(30.0, 150.0),
        filter_order=2,
        zerophase=False,
    ):
        """
        Loads in greens functions and creates all matrices needed

        Args:
            weights: if None, no weighting is applied, array of floats
                corresponding to length and order of st applies manual weighting,
                if 'prenoise' is specified, will used std of noise window before
                event (length by weightpre) or 'distance' to weight by 1/distance
            weightpre (float): length of pre-noise window in seconds (if not None, noise will be used to
                  determine weights)
            period_range (list or tuple): Range of periods to consider in inversion, in seconds
            filter_order (int): Order of filter applied over period_range
            zerophase (bool): If True, zero-phase filtering will be used
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
            self.weight_method = None
        elif isinstance(weights, str):
            # The user specified a weight method
            self.weight_method = weights
        else:
            # The user specified a vector of station weights
            self.weight_method = 'Manual'
            self.weights = weights

        if self.weight_method != 'Manual':
            if weights == 'prenoise' and weightpre is None:
                raise ValueError(
                    'weightpre must be defined if prenoise weighting is used.'
                )
            else:
                self.weightpre = weightpre

        # check if sampling rate specified is compatible with period_range
        if 2.0 * self.filter['freqmax'] > self.sampling_rate:
            raise ValueError(
                'sampling_rate and period_range are not compatible (violates Nyquist).'
            )

        # Always work on copy of data
        st = self.st.copy()

        # filter data to band specified
        st.filter(
            'bandpass',
            freqmin=self.filter['freqmin'],
            freqmax=self.filter['freqmax'],
            corners=self.filter['order'],
            zerophase=self.filter['zerophase'],
        )

        # resample st to sampling_rate
        st.resample(self.sampling_rate)

        # make sure st data are all the same length
        lens = [len(trace.data) for trace in st]
        if len(set(lens)) != 1:
            print(
                'Resampled records are of differing lengths, interpolating all records to same start time and sampling rate'
            )
            stts = [tr.stats.starttime for tr in st]
            lens = [tr.stats.npts for tr in st]
            st.interpolate(
                self.sampling_rate, starttime=np.max(stts), npts=np.min(lens) - 1
            )

        K = 1.0e-15  # CPS variable needed for conversion to meaningful units
        self.datalength = len(st[0].data)

        if self.greenlength > self.datalength:
            raise ValueError(
                'greenlength is greater than datalength. Reselect data and/or '
                'recompute Green\'s functions so that data is longer than Green\'s '
                'functions'
            )

        if self.domain == 'time':
            # ADD WAY TO ACCOUNT FOR WHEN GREENLENGTH IS LONGER THAN DATALENGTH - ACTUALLY SHOULD BE AS LONG AS BOTH ADDED TOGETHER TO AVOID WRAPPING ERROR
            self.lenUall = self.datalength * len(st)
        elif self.domain == 'freq':
            self.NFFT = next_pow_2(
                self.datalength
            )  # +greenlength) #needs to be the length of the two added together because convolution length M+N-1
            self.lenUall = self.NFFT * len(st)
        else:
            raise ValueError('domain not recognized. Must be \'time\' or \'freq\'.')

        if self.method == 'tik':

            self.force_sampling_rate = self.sampling_rate

            # initialize weighting matrices
            Wvec = np.ones(self.lenUall)
            indx = 0
            weight = np.ones(self.numsta)

            n = self.datalength

            for i, trace in enumerate(st):
                # find component of st
                component = trace.stats.channel[2]
                station = trace.stats.station
                if component == 'Z':
                    zvf = read(os.path.join(self.sacdir, '*_%s_*ZVF.sac' % station))
                    if len(zvf) > 1:
                        raise ValueError(f'Found more than one ZVF GF for {station}.')
                    else:
                        zvf = zvf[0]
                    zhf = read(os.path.join(self.sacdir, '*_%s_*ZHF.sac' % station))
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
                        zvff = np.fft.fft(zvf.data, self.NFFT)
                        zhff = np.fft.fft(zhf.data, self.NFFT)
                        ZVF = np.diag(zvff)
                        ZHF = np.diag(zhff)
                    az = math.radians(trace.stats.azimuth)
                    newline = np.hstack(
                        (K * ZVF, K * ZHF * math.cos(az), K * ZHF * math.sin(az))
                    )
                elif component == 'R':
                    rvf = read(os.path.join(self.sacdir, '*_%s_*RVF.sac' % station))
                    if len(rvf) > 1:
                        raise ValueError(f'Found more than one RVF GF for {station}.')
                    else:
                        rvf = rvf[0]
                    rhf = read(os.path.join(self.sacdir, '*_%s_*RHF.sac' % station))
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
                        rvff = np.fft.fft(rvf.data, self.NFFT)
                        rhff = np.fft.fft(rhf.data, self.NFFT)
                        RVF = np.diag(rvff)
                        RHF = np.diag(rhff)
                    az = math.radians(trace.stats.azimuth)
                    newline = np.hstack(
                        (K * RVF, K * RHF * math.cos(az), K * RHF * math.sin(az))
                    )
                elif component == 'T':
                    thf = read(os.path.join(self.sacdir, '*_%s_*THF.sac' % station))
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
                        thff = np.fft.fft(thf.data, self.NFFT)
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
                    datline = np.fft.fft(trace.data, self.NFFT)

                if i == 0:  # initialize G and d if first station
                    G = newline.copy()
                    d = datline.copy()
                else:  # otherwise build on G and d
                    G = np.vstack((G, newline.copy()))
                    d = np.hstack((d, datline.copy()))
                if weights is not None:
                    if self.weight_method == 'Manual':
                        weight[i] = weights[i]
                    elif weights == 'prenoise':
                        weight[i] = 1.0 / np.std(
                            trace.data[0 : int(weightpre * trace.stats.sampling_rate)]
                        )
                    elif weights == 'distance':
                        weight[i] = trace.stats.distance

                    Wvec[indx : indx + self.datalength] = (
                        Wvec[indx : indx + self.datalength] * weight[i]
                    )
                    indx += self.datalength

        elif self.method == 'triangle':

            # initialize weighting matrices
            Wvec = np.ones(self.lenUall)
            indx = 0
            weight = np.ones(self.numsta)

            n = self.datalength
            fshiftby = int(
                self.L / self.sampling_rate
            )  # Number of samples to shift each triangle by
            Flen = int(
                np.floor(self.datalength / fshiftby)
            )  # Number of shifts, corresponds to length of force time function
            self.force_sampling_rate = 1.0 / fshiftby

            for i, trace in enumerate(st):
                # find component of st
                component = trace.stats.channel[2]
                station = trace.stats.station
                if component == 'Z':
                    zvf = read(os.path.join(self.sacdir, '*%s*ZVF.sac' % station))
                    if len(zvf) > 1:
                        raise ValueError(f'Found more than one ZVF GF for {station}.')
                    else:
                        zvf = zvf[0]
                    zhf = read(os.path.join(self.sacdir, '*%s*ZHF.sac' % station))
                    if len(zhf) > 1:
                        raise ValueError(f'Found more than one ZHF GF for {station}.')
                    else:
                        zhf = zhf[0]
                    # process the same way as st (except shouldn't need to resample)
                    """ Don't need to filter these GFs? Has non-zero offset so filtering does weird things, already convolved with LP source-time function
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
                    rvf = read(os.path.join(self.sacdir, '*%s*RVF.sac' % station))
                    if len(rvf) > 1:
                        raise ValueError(f'Found more than one RVF GF for {station}.')
                    else:
                        rvf = rvf[0]
                    rhf = read(os.path.join(self.sacdir, '*%s*RHF.sac' % station))
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
                    thf = read(os.path.join(self.sacdir, '*%s*THF.sac' % station))
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
                    if self.weight_method == 'Manual':
                        weight[i] = weights[i]
                    elif weights == 'prenoise':
                        weight[i] = 1.0 / np.std(
                            trace.data[0 : int(weightpre * trace.stats.sampling_rate)]
                        )
                    elif weights == 'distance':
                        weight[i] = trace.stats.distance

                    Wvec[indx : indx + self.datalength] = (
                        Wvec[indx : indx + self.datalength] * weight[i]
                    )
                    indx += self.datalength

        else:
            raise ValueError(f'Method {self.method} not supported.')

        # Normalize Wvec so largest weight is 1.
        self.Wvec = Wvec / np.max(np.abs(Wvec))
        self.weights = weight / np.max(np.abs(weight))

        if np.shape(G)[0] != len(d):
            raise ValueError(
                'G and d sizes are not compatible, fix something somewhere.'
            )
        self.G = (
            G * 1.0 / self.sampling_rate
        )  # need to multiply G by sample interval (sec) since convolution is an integral
        self.d = d * 100.0  # WHY?convert data from m to cm
        if weights is not None:
            self.W = np.diag(self.Wvec)
        else:
            self.W = None

    def invert(
        self,
        zero_time=None,
        impose_zero=False,
        add_to_zero=False,
        maxduration=None,
        jackknife=False,
        num_iter=200,
        frac_delete=0.5,
        **kwargs,
    ):
        """
        Perform single force inversion of long-period landslide
        seismic signal using Tikhonov regularization

        Args:
            zero_time (float): Optional estimated start time of real part of signal, in seconds from
                start time of seismic data. Useful for making figures showing selected start time
                and also for impose_zero option
            impose_zero (bool): Will add weighting matrix to suggest forces tend towards zero prior
                to zero_time (zero_time must be defined)
            add_to_zero (bool): Add weighting matrix to suggest that all components of force integrate
                to zero.
            maxduration (float): Maximum duration allowed for the event, starting at zero_time if defined,
                otherwise starting from beginning of seismic data. Points after this will tend towards
                zero. This helps tamp down artifacts due to edge effects.
            jackknife
            num_iter
            frac_delete
            kwargs

        Returns: Populates object with the following attributes
            model (array): model vector of concatated components (n x 1) of solution using
                regularization parameter alpha
            Zforce (array): vertical force time series extracted from model
            Nforce (array): same as above for north force
            Eforce (array): same as above for east force
            tvec (array): Time vector, referenced using zero_time (if specified) and corrected for T0
                time shift
            VR (float): Variance reduction (%), rule of thumb, this should be ~50%-80%, if 100%,
                solution is fitting data exactly and results are suspect. If ~5%, model may be wrong or
                something else may be wrong with setup.
            dtorig (array): original data vector (time domain)
            dtnew (array): modeled data vector (Gm-d) (converted to time domain if domain='freq')
            alpha (float): regularization parameter that was used
            fit1 (array):
            size1 (array):
        """

        # Check inputs for consistency
        if impose_zero and not zero_time:
            raise ValueError('impose_zero set to True but no zero_time provided.')

        # Save input choices
        self.regr_param = kwargs  # regression parameters specific to method
        self.add_to_zero = add_to_zero
        self.zero_time = zero_time
        self.impose_zero = impose_zero
        self.maxduration = maxduration

        # Initialize stuff (also serves to clear any previous results if this is a rerun)
        self.model = None
        self.Zforce = None
        self.Nforce = None
        self.Eforce = None
        self.tvec = None
        self.VR = None  # variance reduction
        self.dtorig = None  #
        self.dtnew = None
        self.alpha = None
        self.alphafit = {'alphas': None, 'fit': None, 'size': None}

        if jackknife:
            self.jackknife = dict(
                ZforceL=[],
                ZforceU=[],
                NforceL=[],
                NforceU=[],
                EforceL=[],
                EforceU=[],
                VR_all=[],
                Zforce_all=[],
                Nforce_all=[],
                Eforce_all=[],
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
        """
        Full waveform inversion using Tikhonov regularization

        Args:
            alphaset (float): Set regularization parameter, if None, will search for best alpha using L-curve method
            zero_scaler (float): Factor by which to divide Gnorm to get scaling factor used for zero constraint.
                The lower the number, teh stronger the constraint, but the higher the risk of high freq.
                oscillations due to a sudden release of the constraint
            zero_taper_length (float): length of taper for zero_scaler, in seconds.
                shorter tapers can result in sharp artifacts, longer is better
            tikhonov_ratios (array): Proportion each regularization method contributes, where values correspond
                to [zeroth, first order, second order]. Must add to 1. Only used if method = 'tikh'
        """

        if np.sum(tikhonov_ratios) != 1.0:
            raise ValueError('Tikhonov ratios must add to 1.')
        self.parameters = {}

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

        dl = self.datalength
        gl = int(n / 3)  # self.datalength

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
                    np.floor(((self.zero_time - self.L) * self.force_sampling_rate))
                )  # Potentially need to adjust for T0 here too?
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

        if self.maxduration is not None:
            if self.zero_time is None:
                zerotime = 0.0
            else:
                zerotime = self.zero_time
            startind = int(
                (zerotime + self.T0 + self.maxduration) * self.force_sampling_rate
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
            print('best alpha is %6.1e' % alpha)
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

        if self.domain == 'freq':
            model, residuals, rank, s = sp.linalg.lstsq(A, x)
            self.model = model.copy()
            div = len(model) / 3
            self.Zforce = -np.real(
                np.fft.ifft(model[0:div]) / 10 ** 5
            )  # convert from dynes to newtons, flip so up is positive
            self.Nforce = np.real(np.fft.ifft(model[div : 2 * div]) / 10 ** 5)
            self.Eforce = np.real(np.fft.ifft(model[2 * div :]) / 10 ** 5)
            # run forward model
            df_new = self.G @ model.T
            # convert d and df_new back to time domain
            dt, dtnew = _back2time(self.d, df_new, self.numsta, dl)
            self.dtorig = dt
            self.dtnew = dtnew

        else:  # domain is time
            model, residuals, rank, s = sp.linalg.lstsq(A, x)
            self.model = model.copy()
            div = int(len(model) / 3)
            self.Zforce = (
                -model[0:div] / 10 ** 5
            )  # convert from dynes to netwons, flip so up is positive
            self.Nforce = model[div : 2 * div] / 10 ** 5
            self.Eforce = model[2 * div :] / 10 ** 5
            dtnew = self.G.dot(model)
            self.dtnew = np.reshape(dtnew, (self.numsta, dl))
            self.dtorig = np.reshape(self.d, (self.numsta, dl))

        # compute variance reduction
        self.VR = _varred(self.dtorig, self.dtnew)
        print('variance reduction %f percent' % (self.VR,))
        tvec = (
            np.arange(
                0,
                len(self.Zforce) * 1 / self.force_sampling_rate,
                1 / self.force_sampling_rate,
            )
            - self.T0
        )
        if self.zero_time is not None:
            tvec -= self.zero_time
        if (
            self.method == 'triangle'
        ):  # Shift so that peak of triangle function lines up with time of force interval
            tvec += self.L
        self.tvec = tvec
        self.dtvec = np.arange(
            0, self.datalength / self.sampling_rate, 1 / self.sampling_rate
        )
        if self.zero_time is not None:
            self.dtvec -= self.zero_time
        # Use constant alpha parameter (found above, if not previously set) for jackknife iterations
        stasets = []
        if self.jackknife is not None:
            # Start jackknife iterations
            for ii in range(self.jackknife['num_iter']):
                numcut = int(round(self.jackknife['frac_delete'] * self.numsta))
                numkeep = self.numsta - numcut
                indxcut = rnd.sample(list(range(self.numsta)), numcut)
                stasets.append(indxcut)

                obj = [
                    sum(ind)
                    for ind in zip(
                        np.tile(list(range(self.datalength)), len(indxcut)),
                        np.repeat(
                            [x1 * self.datalength for x1 in indxcut], self.datalength
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

                Aj = Apart + self.alpha ** 2 * (
                    tikhonov_ratios[0] * I
                    + tikhonov_ratios[1] * L1part
                    + tikhonov_ratios[2] * L2part
                )  # Combo of all regularization things (if any are zero they won't matter)
                xj = np.squeeze(Ghat1.conj().T @ dhat1)

                if self.domain == 'freq':
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
                self.jackknife['Zforce_all'].append(Zf.copy())
                self.jackknife['Nforce_all'].append(Nf.copy())
                self.jackknife['Eforce_all'].append(Ef.copy())
                self.jackknife['VR_all'].append(VR.copy())

            self.jackknife['ZforceL'] = np.percentile(
                self.jackknife['Zforce_all'], 2.5, axis=0
            )
            self.jackknife['ZforceU'] = np.percentile(
                self.jackknife['Zforce_all'], 97.5, axis=0
            )
            self.jackknife['EforceL'] = np.percentile(
                self.jackknife['Eforce_all'], 2.5, axis=0
            )
            self.jackknife['EforceU'] = np.percentile(
                self.jackknife['Eforce_all'], 97.5, axis=0
            )
            self.jackknife['NforceL'] = np.percentile(
                self.jackknife['Nforce_all'], 2.5, axis=0
            )
            self.jackknife['NforceU'] = np.percentile(
                self.jackknife['Nforce_all'], 97.5, axis=0
            )

            self.jackknife['VR_all'] = np.array(self.jackknife['VR_all'])

            print(
                'Jackknife VR stats: max %2.0f, min %2.0f, median %2.0f'
                % (
                    self.jackknife['VR_all'].max(),
                    self.jackknife['VR_all'].min(),
                    np.median(self.jackknife['VR_all']),
                )
            )

    def plot_fits(self, equal_scale=True, xlim=None):
        """Create a plot showing the model-produced waveform fit to the data.

        Args:
            equal_scale: If `True`, all plots will share the same y-axis scale
            xlim: Tuple of x-axis limits (time relative to zero time, in seconds)

        Returns:
            The figure handle
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

        for i, tr in enumerate(self.st):
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

    def plotinv(
        self,
        subplots=False,
        xlim=None,
        ylim=None,
        sameY=True,
        highf_tr=None,
        hfylabel=None,
        hfshift=0.0,
        tvecshift=0.0,
        jackshowall=False,
        infra_tr=None,
        infra_shift=0,
    ):
        """
        Plot inversion result

        Args:
            subplots (bool): True, make subplots, False, plot all one one plot
            xlim:
            ylim:
            sameY:
            hfshift:
            infra_shift:
            [ZEN]upper = upper limit of uncertainties (None if none)
            [ZEN]lower = ditto for lower limit
            highf_tr: obspy trace with a start time identical to the start time of the data used in the
                inversion (otherwise won't line up right)
            infra_tr:
            hfylabel (str): Label used for high frequency trace, if not
                defined, will use station name
            tvecshift (float): shift time vector manually by this many seconds,
                not usually need, only for display purposes
            jackshowall (bool): if True and jackknife was run, will show all
                individual runs (will change subplots to True)

        Returns
            figure handle

        """

        tvec = self.tvec - tvecshift

        annot_kwargs = dict(xy=(0.99, 0.25), xycoords='axes fraction', ha='right')

        # Find y limits
        if self.jackknife is None:
            if ylim is None:
                ylim1 = (
                    np.amin([self.Zforce.min(), self.Eforce.min(), self.Nforce.min()]),
                    np.amax([self.Zforce.max(), self.Eforce.max(), self.Nforce.max()]),
                )
                ylim = (
                    ylim1[0] + 0.1 * ylim1[0],
                    ylim1[1] + 0.1 * ylim1[1],
                )  # add 10% on each side to make it look nicer
        else:
            Zupper = self.jackknife['ZforceU']
            Nupper = self.jackknife['NforceU']
            Eupper = self.jackknife['EforceU']
            Zlower = self.jackknife['ZforceL']
            Nlower = self.jackknife['NforceL']
            Elower = self.jackknife['EforceL']
            if ylim is None:
                ylim1 = (
                    np.amin(
                        [
                            Zlower.min(),
                            Elower.min(),
                            Nlower.min(),
                            self.Zforce.min(),
                            self.Nforce.min(),
                            self.Eforce.min(),
                        ]
                    ),
                    np.amax(
                        [
                            Zupper.max(),
                            Eupper.max(),
                            Nupper.max(),
                            self.Zforce.max(),
                            self.Nforce.max(),
                            self.Eforce.max(),
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
                fig = plt.figure(figsize=(10, 8))
                ax1 = fig.add_subplot(311)
                ax2 = fig.add_subplot(312)
                ax3 = fig.add_subplot(313)
            else:
                fig = plt.figure(figsize=(10, 10))
                ax1 = fig.add_subplot(411)
                ax2 = fig.add_subplot(412)
                ax3 = fig.add_subplot(413)
                ax4 = fig.add_subplot(414)

                if infra_tr is not None:
                    plt.close(fig)
                    fig = plt.figure(figsize=(10, 12))
                    ax1 = fig.add_subplot(511)
                    ax2 = fig.add_subplot(512)
                    ax3 = fig.add_subplot(513)
                    ax4 = fig.add_subplot(514)
                    ax5 = fig.add_subplot(515)

            ax1.plot(tvec, self.Zforce, 'b', linewidth=1)
            ax1.set_ylabel('Up force (N)')
            ax2.plot(tvec, self.Nforce, 'r', linewidth=1)
            ax2.set_ylabel('North force (N)')
            ax3.plot(tvec, self.Eforce, 'g', linewidth=1)
            ax3.set_ylabel('East force (N)')

            x = np.concatenate((tvec, tvec[::-1]))
            if self.jackknife is not None:
                if jackshowall:
                    for Z, N, E in zip(
                        self.jackknife['Zforce_all'],
                        self.jackknife['Nforce_all'],
                        self.jackknife['Eforce_all'],
                    ):
                        ax1.plot(self.tvec, Z, 'b', alpha=0.2, linewidth=1)
                        ax2.plot(self.tvec, N, 'r', alpha=0.2, linewidth=1)
                        ax3.plot(self.tvec, E, 'g', alpha=0.2, linewidth=1)
                else:
                    y = np.concatenate((Zlower, Zupper[::-1]))
                    poly = plt.Polygon(
                        list(zip(x, y)), facecolor='b', edgecolor='none', alpha=0.2
                    )
                    ax1.add_patch(poly)
                    y = np.concatenate((Nlower, Nupper[::-1]))
                    poly = plt.Polygon(
                        list(zip(x, y)), facecolor='r', edgecolor='none', alpha=0.2
                    )
                    ax2.add_patch(poly)
                    y = np.concatenate((Elower, Eupper[::-1]))
                    poly = plt.Polygon(
                        list(zip(x, y)), facecolor='g', edgecolor='none', alpha=0.2
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
                        '%s (shifted -%1.0f s)' % (infra_tr.id, infra_shift),
                        **annot_kwargs,
                    )
                else:
                    ax5.annotate('%s' % infra_tr.id, **annot_kwargs)

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
            if self.maxduration is not None:
                [
                    axe.axvline(self.maxduration, color='gray', linestyle='solid', lw=3)
                    for axe in axes
                ]
            if hfylabel is not None:
                ax4.set_ylabel(hfylabel)

            if hfshift != 0:
                ax4.annotate(
                    '%s (shifted -%1.0f s)' % (highf_tr.id, hfshift), **annot_kwargs
                )
            else:
                ax4.annotate('%s' % highf_tr.id, **annot_kwargs)

        else:
            if highf_tr is None:
                fig = plt.figure(figsize=(14, 4))
                ax = fig.add_subplot(111)
            else:
                fig = plt.figure(figsize=(14, 7))
                ax = fig.add_subplot(211)
                ax4 = fig.add_subplot(212)
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
                        '%s - shifted -%1.0f s' % (highf_tr.id, hfshift), **annot_kwargs
                    )
                else:
                    ax4.annotate('%s' % highf_tr.id, **annot_kwargs)

            ax.plot(tvec, self.Zforce, 'b', label='Up')
            ax.plot(tvec, self.Nforce, 'r', label='North')
            ax.plot(tvec, self.Eforce, 'g', label='East')

            if self.jackknife is not None:
                x = np.concatenate((tvec, tvec[::-1]))

                y = np.concatenate((Zlower, Zupper[::-1]))
                poly = plt.Polygon(
                    list(zip(x, y)), facecolor='b', edgecolor='none', alpha=0.2
                )
                ax.add_patch(poly)
                y = np.concatenate((Nlower, Nupper[::-1]))
                poly = plt.Polygon(
                    list(zip(x, y)), facecolor='r', edgecolor='none', alpha=0.2
                )
                ax.add_patch(poly)
                y = np.concatenate((Elower, Eupper[::-1]))
                poly = plt.Polygon(
                    list(zip(x, y)), facecolor='g', edgecolor='none', alpha=0.2
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
            if self.maxduration is not None:
                ax.axvline(self.maxduration, color='gray', linestyle='solid', lw=3)

        t0 = self.st[0].stats.starttime
        if self.zero_time:
            t0 += self.zero_time
        plt.xlabel('Time (s) from {}'.format(t0.strftime('%Y-%m-%d %H:%M:%S')))

        fig.tight_layout()
        fig.show()

        return fig

    def plotangmag(self, xlim=None, ylim=None, tvecshift=0.0):
        """
        plot angles and magnitudes of inversion result and append results
        to object for further use

        USAGE plotinv(Zforce,Nforce,Eforce,tvec,T0,zerotime=0.,subplots=False,Zupper=None,Zlower=None,Eupper=None,Elower=None,Nupper=None,Nlower=None):
        INPUTS
        [ZEN]force
        tvec =
        T0 = T0 (time delay) used in Green's functions (usually negative)
        zerotime = designated time for event start
        vline = plot vertical line at t=vline
        [ZEN]upper = upper limit of uncertainties (None if none)
        [ZEN]lower = ditto for lower limit
        OUPUTS
        fig - figure handle
        """

        tvec = self.tvec - tvecshift

        if self.jackknife is None:
            if ylim is None:
                ylim1 = (
                    np.amin([self.Zforce.min(), self.Eforce.min(), self.Nforce.min()]),
                    np.amax([self.Zforce.max(), self.Eforce.max(), self.Nforce.max()]),
                )
                ylim = (
                    ylim1[0] + 0.1 * ylim1[0],
                    ylim1[1] + 0.1 * ylim1[1],
                )  # add 10% on each side to make it look nicer
        else:
            Zupper = self.jackknife['ZforceU']
            Nupper = self.jackknife['NforceU']
            Eupper = self.jackknife['EforceU']
            Zlower = self.jackknife['ZforceL']
            Nlower = self.jackknife['NforceL']
            Elower = self.jackknife['EforceL']

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
        ax.plot(tvec, self.Zforce, 'b', label='Up')
        ax.plot(tvec, self.Nforce, 'r', label='North')
        ax.plot(tvec, self.Eforce, 'g', label='East')

        if self.jackknife is not None:
            x = np.concatenate((tvec, tvec[::-1]))
            y = np.concatenate((Zlower, Zupper[::-1]))
            poly = plt.Polygon(
                list(zip(x, y)), facecolor='b', edgecolor='none', alpha=0.2
            )
            ax.add_patch(poly)
            y = np.concatenate((Nlower, Nupper[::-1]))
            poly = plt.Polygon(
                list(zip(x, y)), facecolor='r', edgecolor='none', alpha=0.2
            )
            ax.add_patch(poly)
            y = np.concatenate((Elower, Eupper[::-1]))
            poly = plt.Polygon(
                list(zip(x, y)), facecolor='g', edgecolor='none', alpha=0.2
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
        Mag = np.linalg.norm(list(zip(self.Zforce, self.Eforce, self.Nforce)), axis=1)
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
            ax1.plot(tvec, MagL, 'r', label='lower')
            ax1.plot(tvec, MagU, 'r', label='upper')
        else:
            MagU = None
            MagL = None
        ax1.set_ylabel('Force (N)')

        # Plot the horizontal azimuth
        ax2 = fig.add_subplot(413)
        tempang = (180 / np.pi) * np.arctan2(
            self.Nforce, self.Eforce
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

        # Plot the vertical angle
        ax3 = fig.add_subplot(414)
        Vang = (180 / np.pi) * np.arctan(
            self.Zforce / np.sqrt(self.Nforce ** 2 + self.Eforce ** 2)
        )
        ax3.plot(tvec, Vang)
        ax3.set_ylabel('Vertical angle (deg)')

        axes = fig.get_axes()
        if xlim:
            axes = fig.get_axes()
            [axe.set_xlim(xlim) for axe in axes]
            [axe.grid(True) for axe in axes]

        if self.impose_zero:
            [axe.axvline(0, color='gray', linestyle='solid', lw=3) for axe in axes]
        if self.maxduration is not None:
            [
                axe.axvline(self.maxduration, color='gray', linestyle='solid', lw=3)
                for axe in axes
            ]

        plt.xlabel('Time (s)')
        plt.show()

        self.angmag = dict(Mag=Mag, MagU=MagU, MagL=MagL, Vang=Vang, Haz=Haz)

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
        """
        Args:
            filepath (str): full filepath where all files should be saved
                if None, will use self.main_folder
            timestamp (bool): will stamp results with current time so as not
                to overwrite previous results
            figs2save (list): list of figure handles to save
            figs2save_names (list): list of names of figures (appends to end)
            light (bool): to reduce size, does not save seismic data with object
            filetype:

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
            filename = '%s_%1.0f-%1.0fsec_%s%s' % (
                self.nickname,
                self.filter['periodmin'],
                self.filter['periodmax'],
                jk,
                UTCDateTime.now().strftime('%Y-%m-%dT%H%M'),
            )
        else:
            filename = '%s_%1.0f-%1.0fsec_%s' % (
                self.nickname,
                self.filter['periodmin'],
                self.filter['periodmax'],
                jk,
            )

        with open(os.path.join(filepath, '%s.pickle' % filename), 'wb') as f:
            pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

        if figs2save is not None:
            if figs2save_names is None:
                figs2save_names = range(len(figs2save))
            for i, fig in enumerate(figs2save):
                fig.savefig(
                    os.path.join(
                        filepath, '%s_%s.%s' % (filename, figs2save_names[i], filetype)
                    ),
                    bbox_inches='tight',
                )


def _find_alpha(
    Ghat,
    dhat,
    I,
    L1=0.0,
    L2=0.0,
    invmethod='lsq',
    tikhonov_ratios=(1.0, 0.0, 0.0),
    rough=False,
):
    """
    Find best regularization (trade-off) parameter, alpha, by computing model with many values of
    alpha, plotting L-curve, and finding point of steepest curvature where slope is negative.

    Args:
        Ghat (array): m x n matrix of
        dhat (array): 1 x n array of weighted data
        I (array): Identity matrix
        L1 (array): First order roughening matrix, if 0., will use only zeroth order Tikhonov reg.
        L2 (array): Second order roughening matrix, if 0., will use only zeroth order Tikhonov reg.
        invmethod (str): if 'lsq' will use least squares (regular tikhonov), 'nnls' will use non-negative
            least squares
        tikhonov_ratios (list): Proportion each regularization method contributes, where values correspond
            to [zeroth, first order, second order]. Must add to 1.
        rough (bool): If False (default), will do two iterations to fine tune the alpha parameter,
            if True, time will be saved because it will only do one round of searching
    Returns:

    """

    templ1 = np.ceil(
        np.log10(np.linalg.norm(Ghat))
    )  # Roughly estimate largest singular value (should not use alpha larger than expected largest singular value)
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
    # Zero out any points where function is concave so avoid picking points form dropoff at end
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
        # Zero out any points where function is concave so avoid picking points form dropoff at end
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
    """
    Plot L-curve
    """

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(111)
    ax.loglog(fit1, size1, '.')
    for i, alpha in enumerate(alphas):
        text1 = '%3.1E' % alpha
        ax.text(fit1[i], size1[i], text1)
    ax.set_xlabel('Residual norm ||Gm-d||2')
    ax.set_ylabel('Solution norm ||m||2')
    plt.show()
    plt.draw()


def _varred(dt, dtnew):
    """
    compute variance reduction in time domain (%)
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
    """
    convert data back to the time domain and cut off zero padding
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
    """
    Build matrix that can be used for convolution as implemented by matrix multiplication
    size is optional input for desired size as (rows,cols), this will just shift cflip until it
    reaches the right size
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
    """
    Build matrix that can be used for shifting of overlapping triangles for
    triangle method, signal goes across rows and each shift is a new column
    (opposite orientation to _makeconvmat)

    Args:
        c (array): vector of data (usually greens function)
        shiftby (int): number of samples to shift greens function in each row
        size1 (tup): (nrows, ncols) of desired result. Will pad c if nrows is
            greater than len(c). Will shift c forward by shiftby ncols times

    Returns:
        Matrix of shifted c of size size1
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
    """
    Estimate the radius of curvature for each point on line to find corner of L-curve

    Args:
        x (array): x points
        y (array): y points

    Returns:
        radius of curvature for each point (ends will be nan)
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
