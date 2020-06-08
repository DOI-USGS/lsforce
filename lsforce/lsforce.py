import numpy as np
from obspy import Trace, Stream, read, UTCDateTime
from obspy.signal.util import next_pow_2
import math
import glob
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib import cm
import scipy as sp
import urllib.request
import random as rnd
import pickle
from sklearn import linear_model as lm
import xarray as xr
import warnings
import os
import stat
import shutil
import subprocess
import copy
import cartopy.crs as ccrs


class LSForce:
    """
    Class for single force inversions
    """

    def __init__(
        self,
        st,
        samplerate,
        domain='time',
        nickname=None,
        mainfolder=None,
        source_lat=None,
        source_lon=None,
        method='tik',
    ):
        """
        Args:
            st (Stream): obspy stream with event data , source to station
                distance and azimuth must be attached in stats as tr.stats.rdist [km],
                tr.stats.back_azimuth, tr.stats.azimuth. Should be corrected for station
                response but otherwise unfiltered
            samplerate (float): Number of samples per second (Hz) to use in inversion.
                All data will be resampled to this rate and Greens functions
                will be created with this sample rate.
            domain (str): domain in which to do inversion, 'time' (default) or
                'freq'
            nickname (str): Nickname for this event, used for convenient in
                naming files
            mainfolder (str): if None, will use current folder
            source_lat (float): Latitude in decimal degrees of centroid of
                landslide location
            source_lon (float): Longitude in decimal degrees of centroid of
                landslide location
            method (str): 'tik' = full waveform inversion using Tikhonov
                                regularization (L2 norm minimization)
                          'lasso' = full waveform inversion using Lasso method
                                (L1 norm minimization with smoothing)
                          'triangle' = parameterized inversion using overlapping
                                  triangles (variation of method of Ekstrom et al., 2013)
                          'basis' = parameterized using many hanning basis functions
                          'sinusoid' = parameterized using single sinusoid
                                  (variation of method by Chao et al. YEAR)
        """

        # General
        self.st = st
        self.domain = domain
        self.samplerate = samplerate
        self.nickname = nickname
        self.numsta = len(st)
        self.greens_computed = False

        if source_lat is None or source_lon is None:
            raise Exception('source_lat and source_lon not defined')
        else:
            self.lat = source_lat
            self.lon = source_lon

        if mainfolder is None:
            self.mainfolder = os.getcwd()
        else:
            self.mainfolder = mainfolder

        if method not in ['tik', 'lasso', 'triangle']:
            raise Exception('%s method not yet implemented.' % method.upper())

        self.method = method
        if self.method in ['triangle'] and self.domain == 'freq':
            raise Exception(
                'triangle method must be done in time domain, '
                'frequency domain not implemented'
            )

    def compute_greens(self, modelfile, gfduration, T0, L=5.0):
        """
        Use CPS to compute the necessary Greens functions relating the source
        location and the seismic stations in st
        Will compute the type of Greens functions appropriate for the method
        defined during class initiation

        Args:
            modelfile (str): Full file path to location of CPS model file
            gfduration (float):
            T0 (float): number of seconds prior to impulse application
            L (float): half width, in seconds, of triangle. Only needed for
                triangle method, default 5 sec. This will also correspond
                to sample interval of resulting force time function because
                triangles overlap by 50%
        """

        self.modelfile = modelfile
        self.T0 = T0

        self.L = L

        if self.nickname is None:
            self.nickname = ''

        self.moddir = os.path.join(
            self.mainfolder,
            ('%s_%s')
            % (self.nickname, os.path.splitext(os.path.basename(modelfile))[0]),
        )
        self.sacodir = os.path.join(self.moddir, 'sacorig_%s' % self.method)
        self.sacdir = os.path.join(self.moddir, 'sacdata_%s' % self.method)

        # Make all the directories
        if not os.path.exists(self.mainfolder):
            os.mkdir(self.mainfolder)
        if not os.path.exists(self.moddir):
            os.mkdir(self.moddir)
        if not os.path.exists(self.sacodir):
            os.mkdir(self.sacodir)
        if not os.path.exists(self.sacdir):
            os.mkdir(self.sacdir)

        # write T0 file
        with open(os.path.join(self.moddir, 'T0.txt'), 'w') as f:
            f.write(('%3.2f') % T0)

        # write L file, if applicable
        if self.method in ['triangle']:
            with open(os.path.join(self.moddir, 'L.txt'), 'w') as f:
                f.write(('%3.2f') % L)

        # Make sure there is only one occurrence of each station in list (ignore channels)
        stacods = np.unique([tr.stats.station for tr in self.st])
        dists = [self.st.select(station=sta)[0].stats.rdist for sta in stacods]

        # write stadistlist.txt
        f = open(os.path.join(self.moddir, 'stadistlist.txt'), 'w')
        for sta, dis in zip(stacods, dists):
            f.write(('%s\t%5.1f\n') % (sta, dis))
        f.close()

        # write dist file in free format
        # figure out how many samples
        samples = next_pow_2(gfduration * self.samplerate)
        f = open(os.path.join(self.moddir, 'dist'), 'w')
        for dis in dists:
            f.write(
                ('%0.1f %0.2f %i %i 0\n')
                % (dis, 1.0 / self.samplerate, samples, self.T0)
            )
        f.close()
        self.greenlength = samples

        # move copy of modelfile to current directory for recordkeeping
        shutil.copy2(modelfile, os.path.join(self.moddir, os.path.basename(modelfile)))

        # write shell script to run Green's functions
        self.shellscript = os.path.join(self.moddir, 'CPScommands.sh')
        with open(self.shellscript, 'w') as f:
            f.write('#!/bin/bash\n')
            f.write('rm %s\n' % os.path.join(self.sacodir, '*.sac'))
            f.write('rm %s\n' % os.path.join(self.sacdir, '*.sac'))
            f.write(
                'hprep96 -HR 0. -HS 0. -M %s -d %s -R -EXF\n' % (self.modelfile, 'dist')
            )
            f.write('hspec96 > hspec96.out\n')
            if self.method == 'triangle':
                f.write(
                    'hpulse96 -d %s -V -D -t -l %d > Green\n'
                    % ('dist', int(self.L / self.samplerate))
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
            raise Exception('Greens functions were not computed: %s' % stderr)

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
            newname = ('GF_%s_%s_%1.0fkm_%s.sac') % (
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

    def load_greens(self, modelfile):

        """
        If Greens functions were already computed for this exact data
        selection, this simply loads info about them needed for setup

        Args:
            modelfile (str): the name of the model file used to compute the
                Greens functions. This is so they can be found because they
                are saved in a folder referencing the model file name

            TODO add error catching in case the stations in st don't line
            up with computed GFs in folder
        """

        if self.nickname is None:
            self.nickname = ''

        self.moddir = os.path.join(
            self.mainfolder,
            ('%s_%s')
            % (self.nickname, os.path.splitext(os.path.basename(modelfile))[0]),
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
        zeroPhase=False,
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
            zeroPhase (bool): If True, zeroPhase filtering will be used
        """

        # Create filter dictionary to keep track of filter used without
        # creating too many new attributes
        self.filter = {
            'freqmin': 1.0 / period_range[1],
            'freqmax': 1.0 / period_range[0],
            'zeroPhase': zeroPhase,
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
                raise Exception(
                    'weightpre must be defined if prenoise weighting is used'
                )
            else:
                self.weightpre = weightpre

        # check if sampling rate specified is compatible with period_range
        if 2.0 * self.filter['freqmax'] > self.samplerate:
            raise Exception(
                'samplerate and period_range are not compatible, ' 'violates Nyquist'
            )

        # Always work on copy of data
        st = self.st.copy()

        # filter data to band specified
        st.filter(
            'bandpass',
            freqmin=self.filter['freqmin'],
            freqmax=self.filter['freqmax'],
            corners=self.filter['order'],
            zerophase=self.filter['zeroPhase'],
        )

        # resample st to samplerate
        st.resample(self.samplerate)

        # make sure st data are all the same length
        lens = [len(trace.data) for trace in st]
        if len(set(lens)) != 1:
            print(
                'Resampled records are of differing lengths, interpolating all records to same start time and sampling rate'
            )
            stts = [tr.stats.starttime for tr in st]
            lens = [tr.stats.npts for tr in st]
            st.interpolate(
                self.samplerate, starttime=np.max(stts), npts=np.min(lens) - 1
            )

        K = 1.0e-15  # CPS variable needed for conversion to meaningful units
        self.datalength = len(st[0].data)

        if self.greenlength > self.datalength:
            raise Exception(
                'greenlength is greater than datalength. Reselect '
                'data and/or recompute greens functions so that data '
                'is longer than greens functions'
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
            raise Exception('domain not recognized. Must be time or freq')

        if self.method in ['tik', 'lasso']:

            self.Fsamplerate = self.samplerate

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
                        raise Exception('Found more than one ZVF GF for %s' % station)
                    else:
                        zvf = zvf[0]
                    zhf = read(os.path.join(self.sacdir, '*_%s_*ZHF.sac' % station))
                    if len(zhf) > 1:
                        raise Exception('Found more than one ZHF GF for %s' % station)
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
                        zerophase=self.filter['zeroPhase'],
                    )

                    zhf.detrend()
                    zhf.taper(max_percentage=0.05)
                    zhf.filter(
                        'bandpass',
                        freqmin=self.filter['freqmin'],
                        freqmax=self.filter['freqmax'],
                        corners=self.filter['order'],
                        zerophase=self.filter['zeroPhase'],
                    )

                    if self.domain == 'time':
                        ZVF = makeconvmat(zvf.data, size=(n, n))
                        ZHF = makeconvmat(zhf.data, size=(n, n))
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
                        raise Exception('Found more than one RVF GF for %s' % station)
                    else:
                        rvf = rvf[0]
                    rhf = read(os.path.join(self.sacdir, '*_%s_*RHF.sac' % station))
                    if len(rhf) > 1:
                        raise Exception('Found more than one RHF GF for %s' % station)
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
                        zerophase=self.filter['zeroPhase'],
                    )

                    rhf.detrend()
                    rhf.taper(max_percentage=0.05)
                    rhf.filter(
                        'bandpass',
                        freqmin=self.filter['freqmin'],
                        freqmax=self.filter['freqmax'],
                        corners=self.filter['order'],
                        zerophase=self.filter['zeroPhase'],
                    )
                    if self.domain == 'time':
                        RVF = makeconvmat(rvf.data, size=(n, n))
                        RHF = makeconvmat(rhf.data, size=(n, n))
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
                        raise Exception('Found more than one THF GF for %s' % station)
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
                        zerophase=self.filter['zeroPhase'],
                    )
                    if self.domain == 'time':
                        THF = makeconvmat(thf.data, size=(n, n))
                    else:
                        thff = np.fft.fft(thf.data, self.NFFT)
                        THF = np.diag(thff)
                    TVF = 0.0 * THF.copy()
                    az = math.radians(trace.stats.azimuth)
                    newline = np.hstack(
                        (TVF, K * THF * math.sin(az), -K * THF * math.cos(az))
                    )
                else:
                    raise Exception('st not rotated to T and R for %s' % station)
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
                            trace.data[0:int(weightpre * trace.stats.sampling_rate)]
                        )
                    elif weights == 'distance':
                        weight[i] = trace.stats.rdist

                    Wvec[indx:indx + self.datalength] = (
                        Wvec[indx:indx + self.datalength] * weight[i]
                    )
                    indx += self.datalength

        elif self.method in ['triangle']:

            # initialize weighting matrices
            Wvec = np.ones(self.lenUall)
            indx = 0
            weight = np.ones(self.numsta)

            n = self.datalength
            fshiftby = int(
                self.L / self.samplerate
            )  # Number of samples to shift each triangle by
            Flen = int(
                np.floor(self.datalength / fshiftby)
            )  # Number of shifts, corresponds to length of force time function
            self.Fsamplerate = 1.0 / fshiftby

            for i, trace in enumerate(st):
                # find component of st
                component = trace.stats.channel[2]
                station = trace.stats.station
                if component == 'Z':
                    zvf = read(os.path.join(self.sacdir, '*%s*ZVF.sac' % station))
                    if len(zvf) > 1:
                        raise Exception('Found more than one ZVF GF for %s' % station)
                    else:
                        zvf = zvf[0]
                    zhf = read(os.path.join(self.sacdir, '*%s*ZHF.sac' % station))
                    if len(zhf) > 1:
                        raise Exception('Found more than one ZHF GF for %s' % station)
                    else:
                        zhf = zhf[0]
                    # process the same way as st (except shouldn't need to resample)
                    """ Don't need to filter these GFs? Has non-zero offset so filtering does weird things, already convolved with LP source-time function
                    zvf.detrend()
                    zvf.taper(max_percentage=0.05)
                    zvf.filter('bandpass', freqmin=self.filter['freqmin'],
                               freqmax=self.filter['freqmax'],
                               corners=self.filter['order'],
                               zerophase=self.filter['zeroPhase'])

                    zhf.detrend()
                    zhf.taper(max_percentage=0.05)
                    zhf.filter('bandpass', freqmin=self.filter['freqmin'],
                               freqmax=self.filter['freqmax'],
                               corners=self.filter['order'],
                               zerophase=self.filter['zeroPhase'])
                    """
                    ZVF = makeshiftmat(zvf.data, shiftby=fshiftby, size1=(n, Flen))
                    ZHF = makeshiftmat(zhf.data, shiftby=fshiftby, size1=(n, Flen))
                    az = math.radians(trace.stats.azimuth)
                    newline = np.hstack(
                        (K * ZVF, K * ZHF * math.cos(az), K * ZHF * math.sin(az))
                    )
                elif component == 'R':
                    rvf = read(os.path.join(self.sacdir, '*%s*RVF.sac' % station))
                    if len(rvf) > 1:
                        raise Exception('Found more than one RVF GF for %s' % station)
                    else:
                        rvf = rvf[0]
                    rhf = read(os.path.join(self.sacdir, '*%s*RHF.sac' % station))
                    if len(rhf) > 1:
                        raise Exception('Found more than one RHF GF for %s' % station)
                    else:
                        rhf = rhf[0]
                    """ Don't need to filter these GFs?
                    #process the same way as st
                    rvf.detrend()
                    rvf.taper(max_percentage=0.05)

                    rvf.filter('bandpass', freqmin=self.filter['freqmin'],
                               freqmax=self.filter['freqmax'],
                               corners=self.filter['order'],
                               zerophase=self.filter['zeroPhase'])

                    rhf.detrend()
                    rhf.taper(max_percentage=0.05)
                    rhf.filter('bandpass', freqmin=self.filter['freqmin'],
                               freqmax=self.filter['freqmax'],
                               corners=self.filter['order'],
                               zerophase=self.filter['zeroPhase'])
                    """
                    RVF = makeshiftmat(rvf.data, shiftby=fshiftby, size1=(n, Flen))
                    RHF = makeshiftmat(rhf.data, shiftby=fshiftby, size1=(n, Flen))
                    az = math.radians(trace.stats.azimuth)
                    newline = np.hstack(
                        (K * RVF, K * RHF * math.cos(az), K * RHF * math.sin(az))
                    )
                elif component == 'T':
                    thf = read(os.path.join(self.sacdir, '*%s*THF.sac' % station))
                    if len(thf) > 1:
                        raise Exception('Found more than one THF GF for %s' % station)
                    else:
                        thf = thf[0]
                    """ Don't need to filter these GFs?
                    #process the same way as st
                    thf.detrend()
                    thf.taper(max_percentage=0.05)
                    thf.filter('bandpass', freqmin=self.filter['freqmin'],
                               freqmax=self.filter['freqmax'],
                               corners=self.filter['order'],
                               zerophase=self.filter['zeroPhase'])
                    """
                    THF = makeshiftmat(thf.data, shiftby=fshiftby, size1=(n, Flen))
                    TVF = 0.0 * THF.copy()
                    az = math.radians(trace.stats.azimuth)
                    newline = np.hstack(
                        (TVF, K * THF * math.sin(az), -K * THF * math.cos(az))
                    )
                else:
                    raise Exception('st not rotated to T and R for %s' % station)
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
                            trace.data[0: int(weightpre * trace.stats.sampling_rate)]
                        )
                    elif weights == 'distance':
                        weight[i] = trace.stats.rdist

                    Wvec[indx:indx + self.datalength] = (
                        Wvec[indx:indx + self.datalength] * weight[i]
                    )
                    indx += self.datalength

        else:
            # TODO setup for other methods
            print('Put setup for other methods here')

        # Normalize Wvec so largest weight is 1.
        self.Wvec = Wvec / np.max(np.abs(Wvec))
        self.weights = weight / np.max(np.abs(weight))

        if np.shape(G)[0] != len(d):
            raise Exception('G and d sizes are not compatible, fix something somewhere')
        self.G = (
            G * 1.0 / self.samplerate
        )  # need to multiply G by sample interval (sec) since convolution is an integral
        self.d = d * 100.0  # WHY?convert data from m to cm
        if weights is not None:
            self.W = np.diag(self.Wvec)
        else:
            self.W = None

    def invert(
        self,
        zeroTime=None,
        imposeZero=False,
        addtoZero=False,
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
            zeroTime (float): Optional estimated start time of real part of signal, in seconds from
                start time of seismic data. Useful for making figures showing selected start time
                and also for imposeZero option
            imposeZero (bool): Will add weighting matrix to suggest forces tend towards zero prior
                to zeroTime (zeroTime must be defined)
            addtoZero (bool): Add weighting matrix to suggest that all components of force integrate
                to zero.
            maxduration (float): Maximum duration allowed for the event, starting at zeroTime if defined,
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
            tvec (array): Time vector, referenced using zeroTime (if specified) and corrected for T0
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
        if imposeZero and (zeroTime is None or zeroTime == 0.0):
            raise Exception('imposeZero set to True but no zeroTime provided')

        # Save input choices
        self.regr_param = kwargs  # regression parameters specific to method
        self.addtoZero = addtoZero
        self.zeroTime = zeroTime
        self.imposeZero = imposeZero
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

        if self.method in ['tik', 'triangle']:
            self.Tikinvert(**kwargs)
        elif self.method == 'lasso':
            self.Lasso(**kwargs)

    def Tikinvert(
        self,
        alphaset=None,
        alpha_method='Lcurve',
        zeroScaler=15.0,
        zeroTaperlen=20.0,
        Tikhratio=(1.0, 0.0, 0.0),
    ):
        """
        Full waveform inversion using Tikhonov regularization

        Args:
            alphaset (float): Set regularization parameter, if None, will search for best alpha

            alpha_method (str): Method used to find best regularization parameter (alpha) if not defined.
                'Lcurve' chooses based on steepest part of curve and 'Discrepancy' choose based on
                discrepancy principle and noise calculated from data from before zeroTime.
            zeroScaler (float): Factor by which to divide Gnorm to get scaling factor used for zero constraint.
                The lower the number, teh stronger the constraint, but the higher the risk of high freq.
                oscillations due to a sudden release of the constraint
            zeroTaperlen (float): length of taper for zeroScaler, in seconds.
                shorter tapers can result in sharp artifacts, longer is better
            Tikhratio (array): Proportion each regularization method contributes, where values correspond
                to [zeroth, first order, second order]. Must add to 1. Only used if method = 'tikh'
        """

        if np.sum(Tikhratio) != 1.0:
            raise Exception('Tikhonov ratios must add to 1')
        self.parameters = {}

        if self.W is not None:
            Ghat = self.W.dot(self.G)  # np.dot(W.tocsr(),G.tocsr())
            dhat = self.W.dot(self.d)
        else:
            Ghat = self.G  # G.tocsr()
            dhat = self.d

        if self.jackknife is not None:  # save version at this point for use later
            Ghatori = Ghat.copy()
            dhatori = dhat.copy()

        m, n = np.shape(Ghat)

        Ghatnorm = np.linalg.norm(Ghat)

        dl = self.datalength
        gl = int(n / 3)  # self.datalength

        if self.addtoZero is True:  # constrain forces to add to zero
            scaler = Ghatnorm
            first1 = np.hstack((np.ones(gl), np.zeros(2 * gl)))
            second1 = np.hstack((np.zeros(gl), np.ones(gl), np.zeros(gl)))
            third1 = np.hstack((np.zeros(2 * gl), np.ones(gl)))
            A1 = np.vstack((first1, second1, third1)) * scaler
            Ghat = np.vstack((Ghat, A1))
            dhat = np.hstack((dhat, np.zeros(3)))
        else:
            A1 = None

        scaler = Ghatnorm / zeroScaler
        if self.imposeZero:  # tell model when there should be no forces
            # TODO get this to work for triangle method (need to change len methods)
            len2 = int(np.floor(((self.zeroTime + self.T0) * self.Fsamplerate)))
            if self.method == 'triangle':
                len2 = int(
                    np.floor(((self.zeroTime - self.L) * self.Fsamplerate))
                )  # Potentially need to adjust for T0 here too?
            if self.method == 'tik':
                len3 = int(zeroTaperlen * self.Fsamplerate)  # make it constant
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
            if self.zeroTime is None:
                zerotime = 0.0
            else:
                zerotime = self.zeroTime
            startind = int((zerotime + self.T0 + self.maxduration) * self.Fsamplerate)
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
        if Tikhratio[1] != 0.0:
            # Build L1 (first order) roughening matrix
            L1 = np.diag(-1 * np.ones(n)) + np.diag(np.ones(n - 1), k=1)
            L1part = np.dot(L1.T, L1)
        else:
            L1part = 0.0
            L1 = 0.0
        if Tikhratio[2] != 0.0:
            # Build L2 (second order) roughening matrix
            L2 = (
                np.diag(np.ones(n))
                + np.diag(-2 * np.ones(n - 1), k=1)
                + np.diag(np.ones(n - 2), k=2)
            )
            L2part = np.dot(L2.T, L2)
        else:
            L2 = 0.0
            L2part = 0.0

        if alphaset is None:
            if alpha_method == 'Lcurve':
                alpha, fit1, size1, alphas = findalpha(
                    Ghat, dhat, I, L1, L2, Tikhratio=Tikhratio, invmethod='lsq'
                )
            else:
                alpha, fit1, size1, alphas = findalphaD(
                    Ghat,
                    dhat,
                    I,
                    self.zeroTime,
                    self.samplerate,
                    self.numsta,
                    dl,
                    L1=L1,
                    L2=L2,
                    Tikhratio=Tikhratio,
                )  # , tolerance=0.5)
            print('best alpha is %6.1e' % alpha)
            self.alpha = alpha
            self.alphafit['alphas'] = alphas
            self.alphafit['fit'] = fit1
            self.alphafit['size'] = size1
        else:
            self.alpha = alpha

        Ghat = np.matrix(Ghat)
        Apart = np.dot(Ghat.H, Ghat)

        A = Apart + alpha ** 2 * (
            Tikhratio[0] * I + Tikhratio[1] * L1part + Tikhratio[2] * L2part
        )  # Combo of all regularization things (if any are zero they won't matter)
        x = np.squeeze(np.asarray(np.dot(Ghat.H, dhat)))

        if self.domain == 'freq':
            model, residuals, rank, s = sp.linalg.lstsq(A, x)
            self.model = model.copy()
            div = len(model) / 3
            self.Zforce = -np.real(
                np.fft.ifft(model[0:div]) / 10 ** 5
            )  # convert from dynes to newtons, flip so up is positive
            self.Nforce = np.real(np.fft.ifft(model[div:2 * div]) / 10 ** 5)
            self.Eforce = np.real(np.fft.ifft(model[2 * div:]) / 10 ** 5)
            # run forward model
            df_new = np.dot(self.G, model.T)  # forward_model(G,model)
            # convert d and df_new back to time domain
            dt, dtnew = back2time(self.d, df_new, self.numsta, dl)
            self.dtorig = dt
            self.dtnew = dtnew

        else:  # domain is time
            model, residuals, rank, s = sp.linalg.lstsq(A, x)
            self.model = model.copy()
            div = int(len(model) / 3)
            self.Zforce = (
                -model[0:div] / 10 ** 5
            )  # convert from dynes to netwons, flip so up is positive
            self.Nforce = model[div:2 * div] / 10 ** 5
            self.Eforce = model[2 * div:] / 10 ** 5
            dtnew = self.G.dot(model)  # forward_model(G,model)
            self.dtnew = np.reshape(dtnew, (self.numsta, dl))
            self.dtorig = np.reshape(self.d, (self.numsta, dl))

        # compute variance reduction
        self.VR = varred(self.dtorig, self.dtnew)
        print(('variance reduction %f percent') % (self.VR,))
        tvec = (
            np.arange(0, len(self.Zforce) * 1 / self.Fsamplerate, 1 / self.Fsamplerate)
            - self.T0
        )
        if self.zeroTime is not None:
            tvec -= self.zeroTime
        if (
            self.method == 'triangle'
        ):  # Shift so that peak of triangle function lines up with time of force interval
            tvec += self.L
        self.tvec = tvec
        self.dtvec = np.arange(
            0, self.datalength / self.samplerate, 1 / self.samplerate
        )
        if self.zeroTime is not None:
            self.dtvec -= self.zeroTime
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
                Ghat1 = np.matrix(Ghat1)
                Apart = np.dot(Ghat1.H, Ghat1)

                Aj = Apart + self.alpha ** 2 * (
                    Tikhratio[0] * I + Tikhratio[1] * L1part + Tikhratio[2] * L2part
                )  # Combo of all regularization things (if any are zero they won't matter)
                xj = np.squeeze(np.asarray(np.dot(Ghat1.H, dhat1)))

                if self.domain == 'freq':
                    model, residuals, rank, s = sp.linalg.lstsq(Aj, xj)
                    div = len(model) / 3
                    Zf = -np.real(
                        np.fft.ifft(model[0:div]) / 10 ** 5
                    )  # convert from dynes to newtons, flip so up is positive
                    Nf = np.real(np.fft.ifft(model[div:2 * div]) / 10 ** 5)
                    Ef = np.real(np.fft.ifft(model[2 * div:]) / 10 ** 5)
                    # run forward model
                    df_new = np.dot(Gtemp, model.T)  # forward_model(G,model)
                    # convert d and df_new back to time domain
                    dt, dtnew = back2time(dtemp, df_new, numkeep, dl)

                else:  # domain is time
                    model, residuals, rank, s = sp.linalg.lstsq(Aj, xj)
                    div = int(len(model) / 3)
                    Zf = (
                        -model[0:div] / 10 ** 5
                    )  # convert from dynes to netwons, flip so up is positive
                    Nf = model[div: 2 * div] / 10 ** 5
                    Ef = model[2 * div:] / 10 ** 5
                    dtnew = Gtemp.dot(model)  # forward_model(G,model)
                    dt = np.reshape(dtemp, (numkeep, dl))

                VR = varred(dt, dtnew)
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

    def Lasso(
        self,
        G,
        d,
        samplerate,
        numsta,
        datlenorig,
        W=None,
        T0=0,
        alpharatio=10.0,
        domain='time',
        alphaset=None,
        zeroTime=None,
        imposeZero=False,
        addtoZero=False,
        alpha_method='Lcurve',
        maxduration=None,
    ):

        """
        NOT YET UPDATED FOR CLASS STRUCTURE, WONT RUN AS IS
        Wrapper function to perform single force inversion of long-period landslide seismic signal
        using scikit learn's Lasso function (L1 norm minimization) with smoothing added on top

        Args:
            G (array): model matrix (m x n)
            d (array): vector of concatenated data (m x 1)
            samplerate (float): samplerate in Hz of seismic data and green's functions (the two must be equal)
            numsta (int): number of channels used
            datlenorig (int): length, in samples, or original data prior to zero padding etc.
            W: optional weighting matrix, same size as G
            T0 (float): Reference time T0 used in Green's function computation (usually 0.)
            alpharatio (float): Alpha for Lasso will be larger than the alpha for smoothing by a factor
                of alphadiv. If None, only Lasso regularization will be done.
            domain (str): specifies whether calculations are in time domain ('time', default) or
                freq domain ('freq')
            alphaset (float): Set regularization parameter, if None, will search for best alpha
                NOTE ABOUT HOW LASSO HANDLES REG FOR L1 and L2
            zeroTime (float): Optional estimated start time of real part of signal, in seconds from
                start time of seismic data. Useful for making figures showing selected start time
                and also for imposeZero option
            imposeZero (bool): Will add weighting matrix to suggest forces tend towards zero prior
                to zeroTime (zeroTime must be defined)
            addtoZero (bool): Add weighting matrix to suggest that all components of force integrate
                to zero.
            alpha_method (str): Method used to find best regularization parameter (alpha) if not defined.
                'Lcurve' chooses based on steepest part of curve and 'Discrepancy' choose based on
                discrepancy principle and noise calculated from data from before zeroTime.
            maxduration (float): Maximum duration allowed for the event, starting at zeroTime if defined,
                otherwise starting from beginning of seismic data. Points after this will tend towards
                zero. This helps tamp down artifacts due to edge effects.

        Returns: (model, Zforce, Nforce, Eforce, tvec, VR, dt, dtnew, alpha, fit1, size1, alphas, curves)
            model (array): model vector of concatated components (n x 1) of solution using
                regularization parameter alpha
            Zforce (array): vertical force time series extracted from model
            Nforce (array): same as above for north force
            Eforce (array): same as above for east force
            tvec (array): Time vector, referenced using zeroTime (if specified) and corrected for T0
                time shift
            VR (float): Variance reduction (%), rule of thumb, this should be ~50%-80%, if 100%,
                solution is fitting data exactly and results are suspect. If ~5%, model may be wrong or
                something else may be wrong with setup.
            dt (array): original data vector
            dtnew (array): modeled data vector (Gm-d)
            alpha (float): regularization parameter that was used
            fit1 (array):
            size1 (array):
            curves (array):

        """

        raise NotImplementedError('Lasso not yet implemented yet in class structure')

        if W is not None:
            Ghat = W.dot(G)  # np.dot(W.tocsr(),G.tocsr())
            dhat = W.dot(d)
        else:
            Ghat = G  # G.tocsr()
            dhat = d

        m, n = np.shape(Ghat)

        Ghatnorm = np.linalg.norm(Ghat)

        if addtoZero is True:  # constrain forces to add to zero
            scaler = Ghatnorm
            first1 = np.hstack((np.ones(datlenorig), np.zeros(2 * datlenorig)))
            second1 = np.hstack(
                (np.zeros(datlenorig), np.ones(datlenorig), np.zeros(datlenorig))
            )
            third1 = np.hstack((np.zeros(2 * datlenorig), np.ones(datlenorig)))
            A = np.vstack((first1, second1, third1)) * scaler
            Ghat = np.vstack((Ghat, A))
            dhat = np.hstack((dhat, np.zeros(3)))

        scaler = Ghatnorm / 15.0
        if imposeZero is True:  # tell model when there should be no forces
            if zeroTime is None:
                raise Exception('imposeZero set to True but no zeroTime provided')
            len2 = int(np.round((zeroTime) * samplerate))
            len3 = int(
                np.round(0.2 * len2)
            )  # 20% taper overlapping into main event by x seconds
            temp = np.hanning(2 * len3)
            temp = temp[len3:]
            vals = np.hstack((np.ones(len2 - len3), temp))
            for i, val in enumerate(vals):
                first1 = np.zeros(3 * datlenorig)
                second1 = first1.copy()
                third1 = first1.copy()
                first1[i] = val
                second1[i + datlenorig] = val
                third1[i + 2 * datlenorig] = val
                if i == 0:
                    A = np.vstack((first1, second1, third1))
                else:
                    A = np.vstack((A, first1, second1, third1))
            A = A * scaler
            Ghat = np.vstack((Ghat, A))
            dhat = np.hstack((dhat, np.zeros(len(vals) * 3)))

        if maxduration is not None:
            if zeroTime is None:
                zeroTime = 0.0
            startind = int((zeroTime + maxduration) * samplerate)
            len2 = int(datlenorig - startind)
            len3 = int(
                np.round(0.2 * len2)
            )  # 20% taper so zero imposition isn't sudden
            temp = np.hanning(2 * len3)
            temp = temp[:len3]
            vals = np.hstack((temp, np.ones(len2 - len3)))
            for i, val in enumerate(vals):
                place = i + startind
                first1 = np.zeros(3 * datlenorig)
                second1 = first1.copy()
                third1 = first1.copy()
                first1[place] = val
                second1[place + datlenorig] = val
                third1[place + 2 * datlenorig] = val
                if i == 0:
                    A = np.vstack((first1, second1, third1))
                else:
                    A = np.vstack((A, first1, second1, third1))
            A = A * scaler
            Ghat = np.vstack((Ghat, A))
            dhat = np.hstack((dhat, np.zeros(len(vals) * 3)))

        if alphaset is not None:
            alpha = alphaset
        dhat = dhat.T

        # Build roughening matrix
        if alpharatio is not None:
            # Build L2 (second order) roughening matrix
            L2 = (
                np.diag(np.ones(n))
                + np.diag(-2 * np.ones(n - 1), k=1)
                + np.diag(np.ones(n - 2), k=2)
            )

        if alphaset is None:
            if alpharatio is None:
                # Use just LassoCV
                lasso = lm.LassoCV(normalize=False)
                lasso.fit(Ghat, dhat)
                alpha = lasso.alpha_
                fit1 = None
                size1 = None
                alphas = None
                print('best alpha is %6.1e' % alpha)
            else:
                # INSERT STUFF HERE TO FIND ALPHA
                raise Exception(
                    'alphaset=None not implemented yet for smoothed Lasso. Must assign alpha'
                )
        else:
            if alpharatio is not None:
                Ghat2 = np.vstack((Ghat, alphaset * L2))
                dhat2 = np.hstack((dhat, np.zeros(np.shape(L2[0]))))
            else:
                Ghat2 = Ghat
                dhat2 = dhat
                alpharatio = 1.0
            lasso = lm.Lasso(alpha=(alphaset * alpharatio) ** 2)
            lasso.fit(Ghat2, dhat2)
            fit1 = None
            size1 = None
            alphas = None

        model = lasso.coef_
        div = int(len(model) / 3)

        if domain == 'freq':
            Zforce = -np.real(
                np.fft.ifft(model[0:div]) / 10 ** 5
            )  # convert from dynes to newtons, flip so up is positive
            Nforce = np.real(np.fft.ifft(model[div:2 * div]) / 10 ** 5)
            Eforce = np.real(np.fft.ifft(model[2 * div:]) / 10 ** 5)
            # run forward model
            df_new = np.dot(G, model.T)  # forward_model(G,model)
            # convert d and df_new back to time domain
            dt, dtnew = back2time(d, df_new, numsta, datlenorig)

        else:  # domain is time
            Zforce = (
                -model[0:div] / 10 ** 5
            )  # convert from dynes to netwons, flip so up is positive
            Nforce = model[div:2 * div] / 10 ** 5
            Eforce = model[2 * div:] / 10 ** 5
            dtnew = G.dot(model)  # forward_model(G,model)
            dtnew = np.reshape(dtnew, (numsta, datlenorig))
            dt = np.reshape(d, (numsta, datlenorig))

        # compute variance reduction
        VR = varred(dt, dtnew)
        print(('variance reduction %2.0f percent') % (VR,))
        tvec = np.arange(0, len(Zforce) * 1 / samplerate, 1 / samplerate) - T0
        if zeroTime is not None:
            tvec = tvec - zeroTime
        return (
            model,
            Zforce,
            Nforce,
            Eforce,
            tvec,
            VR,
            dt,
            dtnew,
            alpha,
            fit1,
            size1,
            alphas,
        )

    def plotdatafit(self):
        """
        plot comparision between real and synthetic data
        """

        fig = plt.figure(figsize=(10, 11))
        offset = 1.0 * np.amax(np.amax(np.absolute(self.dtorig), 0))
        oneline = offset * np.ones((1, self.datalength))
        labels = []
        yticks1 = []
        for i, trace in enumerate(self.st):
            if i == 0:
                addmat = oneline
            else:
                addmat = np.vstack((addmat, oneline * (i + 1)))
            label = f'{trace.stats.network}.{trace.stats.station} ({trace.stats.channel[-1]}) – {trace.stats.rdist:.1f} km'
            labels.append(label)
            yticks1.append(-(i + 1) * offset)

        ax = fig.add_axes([0.25, 0.05, 0.7, 0.9])
        # .T might be flipping the data upside down...
        ax.plot(
            np.tile(self.dtvec, (self.numsta, 1)).T,
            self.dtorig.T - addmat.T,
            'k',
            label='Original',
        )
        ax.plot(
            np.tile(self.dtvec, (self.numsta, 1)).T,
            self.dtnew.T - addmat.T,
            'r',
            label='Model',
        )
        ax.set_xlim((self.dtvec[0], self.dtvec[-1]))
        ax.set_xlabel('Time (sec)')
        ax.set_yticks(yticks1)
        ax.set_yticklabels(labels)
        ax.set_title('Variance Reduction %2.0f%%' % self.VR)
        redline = mlines.Line2D([], [], color='r', label='Model')
        blackline = mlines.Line2D([], [], color='k', label='Data')
        ax.legend(handles=[blackline, redline], loc='upper right')
        plt.show()
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
                fig = plt.figure(figsize=(14, 9))
                ax1 = fig.add_subplot(311)
                ax2 = fig.add_subplot(312)  # ,sharex=ax1)
                ax3 = fig.add_subplot(313)  # ,sharex=ax1)
            else:
                fig = plt.figure(figsize=(14, 12))
                ax1 = fig.add_subplot(411)
                ax2 = fig.add_subplot(412)  # ,sharex=ax1)
                ax3 = fig.add_subplot(413)  # ,sharex=ax1)
                ax4 = fig.add_subplot(414)

                if infra_tr is not None:
                    plt.close(fig)
                    fig = plt.figure(figsize=(14, 15))
                    ax1 = fig.add_subplot(511)
                    ax2 = fig.add_subplot(512)  # ,sharex=ax1)
                    ax3 = fig.add_subplot(513)  # ,sharex=ax1)
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
                    raise Exception('highf_tr is not an obspy trace')
                tvec2 = np.linspace(
                    0,
                    (len(highf_tr.data) - 1) * 1 / highf_tr.stats.sampling_rate,
                    num=len(highf_tr.data),
                )
                # Temporary fix, adjust for same zerotime
                if self.zeroTime:
                    tvec2 -= self.zeroTime
                tvec2 -= hfshift
                ms2ums = 1e6
                ax4.plot(tvec2, highf_tr.data * ms2ums, 'black')
                ax4.set_ylabel('Velocity (μm/s)')

            if infra_tr is not None:
                if type(infra_tr) != Trace:
                    raise Exception('highf_tr is not an obspy trace')
                tvec2 = np.linspace(
                    0,
                    (len(infra_tr.data) - 1) * 1 / infra_tr.stats.sampling_rate,
                    num=len(infra_tr.data),
                )
                # Temporary fix, adjust for same zerotime
                if self.zeroTime:
                    tvec2 -= self.zeroTime
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

            if self.imposeZero:
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
                    raise Exception('highf_tr is not an obspy trace')
                tvec2 = np.linspace(
                    0,
                    (len(highf_tr.data) - 1) * 1 / highf_tr.stats.sampling_rate,
                    num=len(highf_tr.data),
                )
                # Temporary fix, adjust for same zerotime
                if self.zeroTime:
                    tvec2 -= self.zeroTime
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
            if self.imposeZero:
                ax.axvline(0, color='gray', linestyle='solid', lw=3)
            if self.maxduration is not None:
                ax.axvline(self.maxduration, color='gray', linestyle='solid', lw=3)

        t0 = self.st[0].stats.starttime
        if self.zeroTime:
            t0 += self.zeroTime
        plt.xlabel('Time (s) from {}'.format(t0.strftime('%Y-%m-%d %H:%M:%S')))
        plt.show()
        return fig

    def plotangmag(
        self, xlim=None, ylim=None, tvecshift=0.0
    ):
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

        if self.imposeZero:
            [axe.axvline(0, color='gray', linestyle='solid', lw=3) for axe in axes]
        if self.maxduration is not None:
            [
                axe.axvline(self.maxduration, color='gray', linestyle='solid', lw=3)
                for axe in axes
            ]

        plt.xlabel('Time (sec')
        plt.show()

        self.angmag = dict(Mag=Mag, MagU=MagU, MagL=MagL, Vang=Vang, Haz=Haz)

        return fig

    def _integrate_acceleration(
        self, Zforce, Eforce, Nforce, Mass, startidx, endidx, detrend=None
    ):

        traj_tvec = self.tvec[startidx:endidx + 1]

        dx = 1.0 / self.Fsamplerate
        Za = -Zforce.copy()[startidx:endidx + 1] / Mass
        Ea = -Eforce.copy()[startidx:endidx + 1] / Mass
        Na = -Nforce.copy()[startidx:endidx + 1] / Mass
        Zvel = np.cumsum(Za) * dx
        Evel = np.cumsum(Ea) * dx
        Nvel = np.cumsum(Na) * dx

        # Detrend is either None (no detrending) or a time where velo should
        # be fully tapered to zero
        if detrend:
            zeroidx = np.where(traj_tvec == detrend)[0][
                0
            ]  # Index corresponding to time where velo should be zero
            for comp in [Zvel, Evel, Nvel]:
                trend = np.linspace(0, comp[zeroidx], len(traj_tvec[:zeroidx]))
                comp[:zeroidx] -= trend
                comp[zeroidx:] = np.zeros(len(comp[zeroidx:]))

        Zdisp = np.cumsum(Zvel) * dx
        Edisp = np.cumsum(Evel) * dx
        Ndisp = np.cumsum(Nvel) * dx

        return Za, Ea, Na, Zvel, Evel, Nvel, Zdisp, Edisp, Ndisp, traj_tvec

    def _trajectory_automass(
        self,
        Zforce,
        Eforce,
        Nforce,
        Mass=None,
        target_length=None,
        duration=None,
        detrend=None,
    ):
        """
        Calls _integrate_acceleration().
        """

        # Check args
        if Mass and target_length:
            raise ValueError('Cannot specify both mass and target length!')
        if not Mass and not target_length:
            raise ValueError('You must specify either mass or target length!')

        startidx = np.where(self.tvec == 0)[0][0]  # Always start at t = 0
        # Clip time series to `duration` [s] if desired
        if duration:
            endidx = np.where(self.tvec == duration)[0][0]
        else:
            endidx = len(self.tvec)

        # Either use the mass that was provided, or calculate one
        if target_length:

            MASS_INC = int(1e7)  # [kg] Smaller increment is slower but more precise

            # Initialize with end-members
            Mass = 0  # [kg]
            current_length = np.inf  # [km]

            while current_length > target_length:

                Mass += MASS_INC  # Increase the mass

                # Calculate the runout length [km] based on this mass
                *_, Edisp, Ndisp, _ = self._integrate_acceleration(
                    Zforce, Eforce, Nforce, Mass, startidx, endidx, detrend
                )
                current_length = _calculate_Hdist(Edisp, Ndisp)[-1] / 1000  # [km]
        else:
            Mass = int(Mass)

        # Calculate trajectory based on mass assigned above
        (
            Za,
            Ea,
            Na,
            Zvel,
            Evel,
            Nvel,
            Zdisp,
            Edisp,
            Ndisp,
            traj_tvec,
        ) = self._integrate_acceleration(
            Zforce, Eforce, Nforce, Mass, startidx, endidx, detrend
        )

        return Za, Ea, Na, Zvel, Evel, Nvel, Zdisp, Edisp, Ndisp, Mass, traj_tvec

    def trajectory(
        self,
        Mass=None,
        target_length=None,
        duration=None,
        elevation_profile=False,
        plot_jackknife=False,
        image=None,
        dem=None,
        reference_point=None,
        detrend_velocity=None,
    ):
        """
        Integrate force time series to velocity and then displacement. Either
        provide a mass or a target horizontal runout length. If a length is
        provided, the code will find the mass that achieves this length. Calls
        _trajectory_automass().

        Args:
            Mass: Landslide mass [kg]
            target_length: Horizontal runout length from groundtruth [m]
            duration: Clip time series to go from 0-duration [s]
            elevation_profile: If True, plot vertical displacement versus
                               horizontal runout distance (H vs. L) instead of
                               a map view
            plot_jackknife: Toggle plotting jackknifed displacements as well
            image: An xarray.DataArray with coordinates defined in km with the
                   origin (0, 0) being the start location of the trajectory
            dem: A UTM-projected DEM GeoTIFF to slice thru for elevation
                 profile plot
            reference_point (int/float or list): Plot a dot on trajectory, and
                                                 line on colorbar, at this
                                                 specified time(s) for
                                                 reference (default: None, for
                                                 no markings)
            detrend_velocity: If provided, force velocity to linearly go to
                              zero at this time [s]. If None, don't detrend

        Returns:
            The output figure
        """

        # Convert reference points to numpy array
        reference_points = np.atleast_1d(reference_point)

        # For the full inversion (all channels) result
        (
            self.Za,
            self.Ea,
            self.Na,
            self.Zvel,
            self.Evel,
            self.Nvel,
            self.Zdisp,
            self.Edisp,
            self.Ndisp,
            self.Mass,
            self.traj_tvec,
        ) = self._trajectory_automass(
            self.Zforce,
            self.Eforce,
            self.Nforce,
            Mass=Mass,
            target_length=target_length,
            duration=duration,
            detrend=detrend_velocity,
        )
        self.Hdist = _calculate_Hdist(self.Edisp, self.Ndisp)

        fig, ax = plt.subplots()

        # Converting to km below
        if elevation_profile:
            x = self.Hdist / 1000
            y = self.Zdisp / 1000
        else:
            x = self.Edisp / 1000
            y = self.Ndisp / 1000
        sc = ax.scatter(x, y, c=self.traj_tvec, cmap='rainbow', zorder=100)

        if elevation_profile:
            ax.set_xlabel('Horizontal distance (km)')
            ax.set_ylabel('Vertical distance (km)')
        else:
            ax.set_xlabel('East distance (km)')
            ax.set_ylabel('North distance (km)')

        t0 = self.st[0].stats.starttime
        if self.zeroTime:
            t0 += self.zeroTime
        cbar = plt.colorbar(
            sc, label='Time (s) from {}'.format(t0.strftime('%Y-%m-%d %H:%M:%S'))
        )

        # Plot reference points, if any
        if np.any(reference_points):
            cmap = cm.get_cmap('Greys_r', reference_points.size)
            for i, time in enumerate(reference_points):
                try:
                    ref_pt_ind = np.where(self.traj_tvec == time)[0][0]
                except IndexError:
                    raise  # No point corresponding to requested reference time
                ax.scatter(x[ref_pt_ind], y[ref_pt_ind], color=cmap(i), zorder=150)
                cbar.ax.plot(
                    [self.traj_tvec.min(), self.traj_tvec.max()],
                    [time, time],
                    color=cmap(i),
                    linewidth=2,
                )

        title = f'mass = {self.Mass:,} kg\nrunout length = {self.Hdist[-1]/1000:.2f} km'
        if target_length:
            title += f'\n(target length = {target_length:g} km)'
        ax.set_title(title)

        # Plot jackknife trajectories as well if desired
        if plot_jackknife:
            self.jackknife['Zdisp_all'] = []
            self.jackknife['Edisp_all'] = []
            self.jackknife['Ndisp_all'] = []
            self.jackknife['Hdist_all'] = []
            for i in range(self.jackknife['num_iter']):
                *_, Zdisp_i, Edisp_i, Ndisp_i, _, _ = self._trajectory_automass(
                    self.jackknife['Zforce_all'][i],
                    self.jackknife['Eforce_all'][i],
                    self.jackknife['Nforce_all'][i],
                    Mass=Mass,
                    target_length=target_length,
                    duration=duration,
                    detrend=detrend_velocity,
                )

                # Store jackknifed trajectories
                self.jackknife['Zdisp_all'].append(Zdisp_i)
                self.jackknife['Edisp_all'].append(Edisp_i)
                self.jackknife['Ndisp_all'].append(Ndisp_i)

                # Converting to km below
                if elevation_profile:
                    Hdist_i = _calculate_Hdist(Edisp_i, Ndisp_i)
                    self.jackknife['Hdist_all'].append(Hdist_i)  # Store
                    x = Hdist_i / 1000
                    y = Zdisp_i / 1000
                else:
                    x = Edisp_i / 1000
                    y = Ndisp_i / 1000
                ax.scatter(x, y, c=self.traj_tvec, cmap='rainbow', alpha=0.02)

        ax.axis('equal')

        if (image is not None) and (not elevation_profile):
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
            image.plot.imshow(
                ax=ax, cmap='Greys_r', add_colorbar=False, add_labels=False, zorder=-10
            )
            ax.axis('equal')
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)

        if (dem is not None) and elevation_profile:
            distance, drop = self._slice_dem(dem)
            ax.plot(distance / 1000, drop / 1000, color='black', zorder=100)

        plt.tight_layout()
        plt.show()

        return fig

    def _slice_dem(self, dem_file, interp_spacing=0.1):
        """
        Slice through an input DEM along the trajectory path.

        Args:
            dem_file: DEM GeoTIFF to slice (must be UTM-projected!)
            interp_spacing: [m] Density of interpolation points

        Returns:
            horizontal_distance: Distance along path given by EDisp, NDisp
            elevation: Elevation along horizontal_distance
        """

        dem = xr.open_rasterio(dem_file).squeeze()
        # Set no data values to NaN
        dem = dem.where(dem != dem.nodatavals)

        # Define interpolation points in UTM space
        crs = ccrs.epsg(int(dem.crs.split(':')[-1]))
        if crs.proj4_params['proj'] != 'utm':
            raise ValueError('Input DEM must have a UTM projection!')
        loc_utm = crs.transform_point(self.lon, self.lat, ccrs.Geodetic())
        points = [
            [x + loc_utm[0], y + loc_utm[1]] for x, y in zip(self.Edisp, self.Ndisp)
        ]

        # Densify the coarse points
        path_x, path_y = [], []
        x_prev, y_prev = points[0]
        for pt in points[1:]:
            x, y = pt

            seg_length = np.linalg.norm([y - y_prev, x - x_prev])
            # Choose a number of pts that gives approximately interp_spacing
            n = int(seg_length / interp_spacing)

            # Append densified path
            path_x = np.hstack([path_x, np.linspace(x_prev, x, n)])
            path_y = np.hstack([path_y, np.linspace(y_prev, y, n)])

            x_prev, y_prev = x, y

        # Actually interpolate!
        profile = dem.interp(
            x=xr.DataArray(path_x), y=xr.DataArray(path_y), method='linear'
        )

        # Find horizontal distance vector (works for curvy paths!)
        horiz_dist = np.hstack(
            [
                0,
                np.cumsum(
                    np.linalg.norm([np.diff(profile.x), np.diff(profile.y)], axis=0)
                ),
            ]
        )

        # Check that interp_spacing wasn't too coarse by matching path lengths
        if not np.isclose(horiz_dist[-1], self.Hdist[-1]):
            raise ValueError('interp_spacing was too coarse. Try decreasing.')

        warnings.warn('Assuming DEM vertical unit is meters!')

        profile.data -= profile.data[0]  # Start at 0 and go negative

        return horiz_dist, profile.data

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
                if None, will use self.mainfolder
            timestamp (bool): will stamp results with current time so as not
                to overwrite previous results
            figs2save (list): list of figure handles to save
            figs2save_names (list): list of names of figures (appends to end)
            light (bool): to reduce size, does not save seismic data with object
            filetype:

        """
        if filepath is None:
            filepath = self.mainfolder

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


def _calculate_Hdist(Edisp, Ndisp):
    """
    Calculate horizontal distance vector (Hdist) from east and north
    displacement vectors. This is the horizontal distance "along the avalanche
    path" as a function of time. Hdist[-1] is L, the horizontal runout distance
    (which is shorter than the 3-D runout distance).

    Args:
        Edisp: Eastward displacement vector as a function of time [m]
        Ndisp: Northward displacement vector as a function of time [m]

    Returns:
        Horizontal distance as a function of time [m]
    """

    dx = np.diff(Edisp)
    dy = np.diff(Ndisp)

    return np.hstack([0, np.cumsum(np.linalg.norm([dx, dy], axis=0))])


def findalpha(
    Ghat,
    dhat,
    I,
    L1=0.0,
    L2=0.0,
    invmethod='lsq',
    Tikhratio=(1.0, 0.0, 0.0),
    rough=False,
):
    """
    Find best regularization (trade-off) parameter, alpha, by computing model with many values of
    alpha, plotting Lcurve, and finding point of steepest curvature where slope is negative.

    Args:
        Ghat (array): m x n matrix of
        dhat (array): 1 x n array of weighted data
        I (array): Identity matrix
        L1 (array): First order roughening matrix, if 0., will use only zeroth order Tikhonov reg.
        L2 (array): Second order roughening matrix, if 0., will use only zeroth order Tikhonov reg.
        invmethod (str): if 'lsq' will use least squares (regular tikhonov), 'nnls' will use non-negative
            least squares
        Tikhratio (list): Proportion each regularization method contributes, where values correspond
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

    # Convert to matrix so can use complex conjugate .H (so this code will work for both time and freq domains)
    Ghat = np.matrix(Ghat)

    Apart = np.dot(Ghat.H, Ghat)
    if type(L2) == float or type(L2) == int:
        L2part = 0.0
    else:
        L2part = np.dot(L2.T, L2)
    if type(L1) == float or type(L1) == int:
        L1part = 0.0
    else:
        L1part = np.dot(L1.T, L1)

    x = np.squeeze(np.asarray(np.dot(Ghat.H, dhat)))

    # rough first iteration
    for alpha in alphas:
        A = Apart + alpha ** 2 * (
            Tikhratio[0] * I + Tikhratio[1] * L1part + Tikhratio[2] * L2part
        )  # Combo of all regularization things
        if invmethod == 'lsq':
            model, residuals, rank, s = sp.linalg.lstsq(A, x)
        elif invmethod == 'nnls':
            model, residuals = sp.optimize.nnls(A, x)
        else:
            raise Exception('inversion method %s not recognized' % invmethod)
        temp1 = np.dot(Ghat, model.T) - dhat
        fit1.append(sp.linalg.norm(temp1))
        size1.append(
            sp.linalg.norm(Tikhratio[0] * model)
            + sp.linalg.norm(Tikhratio[1] * np.dot(L1part, model))
            + sp.linalg.norm(Tikhratio[2] * np.dot(L2part, model))
        )
    fit1 = np.array(fit1)
    size1 = np.array(size1)
    curves = curvature(np.log10(fit1), np.log10(size1))
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
                Tikhratio[0] * I + Tikhratio[1] * L1part + Tikhratio[2] * L2part
            )  # Combo of all regularization things
            if invmethod == 'lsq':
                model, residuals, rank, s = sp.linalg.lstsq(A, x)
            elif invmethod == 'nnls':
                model, residuals = sp.optimize.nnls(A, x)

            temp1 = np.dot(Ghat, model.T) - dhat
            fit1.append(sp.linalg.norm(temp1))
            size1.append(
                sp.linalg.norm(Tikhratio[0] * model)
                + sp.linalg.norm(Tikhratio[1] * np.dot(L1part, model))
                + sp.linalg.norm(Tikhratio[2] * np.dot(L2part, model))
            )
        fit1 = np.array(fit1)
        size1 = np.array(size1)
        curves = curvature(np.log10(fit1), np.log10(size1))
        # Zero out any points where function is concave so avoid picking points form dropoff at end
        slp2 = np.gradient(np.gradient(np.log10(size1), np.log10(fit1)), np.log10(fit1))
        alphas = np.array(alphas)
        tempcurve = curves.copy()
        tempcurve[slp2 < 0] = np.max(curves)
        idx = np.argmin(tempcurve)
        bestalpha = alphas[idx]
    else:
        bestalpha = alpha

    Lcurve(fit1, size1, alphas)
    if type(bestalpha) == list:
        if len(bestalpha) > 1:
            raise Exception('Returned more than one alpha value, check codes')
        bestalpha = bestalpha[0]
    return bestalpha, fit1, size1, alphas


def findalphaD(
    Ghat,
    dhat,
    I,
    zeroTime,
    samplerate,
    numsta,
    datlenorig,
    tolerance=None,
    L1=0.0,
    L2=0.0,
    Tikhratio=(1.0, 0.0, 0.0),
):
    """
    Use discrepancy principle and noise window to find best alpha (tends to find value that
    is too large such that amplitudes are damped)

    Uses Mozorokov discrepancy principle (and bisection method in log-log space?) to find an
    appropriate value of alpha that results in a solution with a fit that is slightly larger than
    the estimated noise level

    Args:
        tolerance (float): how close you want to get to the noise level with the solution
        Ghat:
        dhat:
        I:
        zeroTime:
        samplerate:
        numsta:
        datlenorig:
        L1:
        L2:
        Tikhratio:
    """

    # Estimate the noise level (use signal before zeroTime)
    dtemp = dhat.copy()[: numsta * datlenorig]  # Trim off any extra zeros
    dtemp = np.reshape(dtemp, (numsta, datlenorig))
    if zeroTime is None:
        print('zeroTime not defined, noise estimated from first 100 samples')
        samps = 100
    else:
        samps = int(zeroTime * samplerate)
    temp = dtemp[:, :samps]
    noise = np.sum(np.sqrt(datlenorig * np.std(temp, axis=1) ** 2))

    # Find ak and bk that yield f(alpha) = ||Gm-d|| - ||noise|| that have f(alpha) values with
    # opposite signs
    templ1 = np.floor(np.log10(np.linalg.norm(Ghat)))
    templ2 = np.arange(templ1 - 6, templ1)
    ak = 10 ** templ2[0]
    bk = 10 ** templ2[-1]
    opposite = False
    fit1 = []
    size1 = []
    alphas = []

    Ghat = np.matrix(Ghat)

    Apart = np.dot(Ghat.H, Ghat)
    if type(L2) == float or type(L2) == int:
        L2part = 0.0
    else:
        L2part = np.dot(L2.T, L2)
    if type(L1) == float or type(L2) == int:
        L1part = 0.0
    else:
        L1part = np.dot(L1.T, L1)

    x = np.squeeze(np.asarray(np.dot(Ghat.H, dhat)))

    while opposite is False:
        print(('ak = %s' % (ak,)))
        print(('bk = %s' % (bk,)))
        Aa = Apart + ak ** 2 * (
            Tikhratio[0] * I + Tikhratio[1] * L1part + Tikhratio[2] * L2part
        )
        modelak, residuals, rank, s = sp.linalg.lstsq(Aa, x)
        Ab = Apart + bk ** 2 * (
            Tikhratio[0] * I + Tikhratio[1] * L1part + Tikhratio[2] * L2part
        )
        modelbk, residuals, rank, s = sp.linalg.lstsq(Ab, x)
        fitak = sp.linalg.norm(np.dot(Ghat, modelak.T) - dhat)
        fitbk = sp.linalg.norm(np.dot(Ghat, modelbk.T) - dhat)
        # Save info on these runs for Lcurve later if desired
        fit1.append(fitak)
        alphas.append(ak)
        size1.append(
            sp.linalg.norm(Tikhratio[0] * modelak)
            + sp.linalg.norm(Tikhratio[1] * np.dot(L1part, modelak))
            + sp.linalg.norm(Tikhratio[2] * np.dot(L2part, modelak))
        )
        fit1.append(fitbk)
        alphas.append(bk)
        size1.append(
            sp.linalg.norm(Tikhratio[0] * modelbk)
            + sp.linalg.norm(Tikhratio[1] * np.dot(L1part, modelbk))
            + sp.linalg.norm(Tikhratio[2] * np.dot(L2part, modelbk))
        )
        fak = fitak - noise  # should be negative
        fbk = fitbk - noise  # should be positive
        print(('fak = %s' % (fak,)))
        print(('fbk = %s' % (fbk,)))
        if fak * fbk < 0:
            opposite = True
        if fak > 0:
            ak = 10 ** (np.log10(ak) - 1)
        if fbk < 0:
            bk = 10 ** (np.log10(bk) + 1)

    if tolerance is None:
        tolerance = noise / 10.0

    # Now use bisection method to find the best alpha value within tolerance
    tol = noise + 100.0
    while tol > tolerance:
        # Figure out whether to change ak or bk
        # Compute midpoint (in log units)
        ck = 10 ** (0.5 * (np.log10(ak) + np.log10(bk)))
        Ac = Apart + ck ** 2 * (
            Tikhratio[0] * I + Tikhratio[1] * L1part + Tikhratio[2] * L2part
        )
        modelck, residuals, rank, s = sp.linalg.lstsq(Ac, x)
        fitck = sp.linalg.norm(np.dot(Ghat, modelck.T) - dhat)
        fit1.append(fitck)
        alphas.append(ck)
        size1.append(
            sp.linalg.norm(Tikhratio[0] * modelck)
            + sp.linalg.norm(Tikhratio[1] * np.dot(L1part, modelck))
            + sp.linalg.norm(Tikhratio[2] * np.dot(L2part, modelck))
        )
        fck = fitck - noise
        print(('ck = %s' % (ck,)))
        print(('fitck = %s' % (fitck,)))
        print(('fck = %s' % (fck,)))
        tol = np.abs(fck)
        if fck * fak < 0:
            bk = ck
        else:
            ak = ck
        print(('ak = %s' % (ak,)))
        print(('bk = %s' % (bk,)))

    bestalpha = ck
    print(('best alpha = %s' % (bestalpha,)))
    fit1 = np.array(fit1)
    size1 = np.array(size1)
    Lcurve(fit1, size1, alphas)
    return bestalpha, fit1, size1, alphas


def Lcurve(fit1, size1, alphas):
    """
    Plot Lcurve
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


def varred(dt, dtnew):
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


def back2time(d, df_new, numsta, datlenorig):
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


def forward_model(G, model):
    """
    run the forward model (without weights in order to compare to unweighted data)
    """

    dnew = np.dot(G, model.T)
    return dnew


def makeconvmat(c, size=None):
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
                p = np.concatenate(((cflip[-(i + 1):]), np.zeros(2 * len(c))))
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
                p = np.concatenate(((cflip[-(i + 1):]), np.zeros(size[1])))
            p = p[: size[1]]  # cut p to the right size
            C[i, :] = p.copy()
    return C


def makeshiftmat(c, shiftby, size1):
    """
    Build matrix that can be used for shifting of overlapping triangles for
    triangle method, signal goes across rows and each shift is a new column
    (opposite orientation to makeconvmat)

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


def curvature(x, y):
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
        xsub = x[i - 1:i + 2]
        ysub = y[i - 1:i + 2]
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


def rotate2ZNE(
    data_1, azimuth_1, dip_1, data_2, azimuth_2, dip_2, data_3, azimuth_3, dip_3
):
    """
    taken from https://github.com/obspy/obspy/blob/a8cf88bfc28b7d06b88427ca626f67418922aa52/obspy/signal/rotate.py#L188-251Rotates an arbitrarily oriented three-component vector to ZNE.
    Each components orientation is described with a azimuth and a dip. The
    azimuth is defined as the degrees from North, clockwise and the dip is the
    defined as the number of degrees, down from horizontal. Both definitions
    are according to the SEED standard.
    The three components need not be orthogonal to each other but the
    components have to be linearly independent. The function performs a full
    base change to orthogonal Vertical, North, and East orientations.
    :param data_1: Data component 1.
    :param azimuth_1: The azimuth of component 1.
    :param dip_1: The dip of component 1.
    :param data_2: Data component 2.
    :param azimuth_2: The azimuth of component 2.
    :param dip_2: The dip of component 2.
    :param data_3: Data component 3.
    :param azimuth_3: The azimuth of component 3.
    :param dip_3: The dip of component 3.
    :rtype: Tuple of three NumPy arrays.
    :returns: The three rotated components, oriented in Z, N, and E.
    """

    # Internally works in Vertical, South, and East components; a right handed
    # coordinate system.
    # Define the base vectors of the old base in terms of the new base vectors.
    base_vector_1 = _dip_azimuth2ZSE_base_vector(dip_1, azimuth_1)
    base_vector_2 = _dip_azimuth2ZSE_base_vector(dip_2, azimuth_2)
    base_vector_3 = _dip_azimuth2ZSE_base_vector(dip_3, azimuth_3)
    # Build transformation matrix.
    T = np.matrix([base_vector_1, base_vector_2, base_vector_3]).transpose()
    # Apply it.
    z, s, e = np.dot(T, [data_1, data_2, data_3])
    # Replace all negative zeros. These might confuse some further processing
    # programs.
    z = np.array(z).ravel()
    z[z == -0.0] = 0
    # Return a North pointing array.
    n = -1.0 * np.array(s).ravel()
    n[n == -0.0] = 0
    e = np.array(e).ravel()
    e[e == -0.0] = 0
    return z, n, e


def rotate(st, baz=None):
    """
    rotate all components of st that can be rotated, to radial and transverse

    Args:
        st: obspy stream object to rotate
        baz (list): Not required if backaz already attached to st stats, list
            of backazimuths corresponding to st
    """

    # implant baz in st's
    if baz:
        for i, trace in enumerate(st):
            trace.stats.back_azimuth = baz[i]
            st[i] = trace
    else:
        try:
            st[0].stats.back_azimuth
        except:
            print('need to attach baz')
            return
    # get list of station location code pairs present
    staloc = list(
        set([trace.stats.station + '.' + trace.stats.location for trace in st])
    )
    st_rotated = Stream(traces=Trace())  # initialize, pop this one off later
    for station in staloc:
        # get all components with that station name
        loc = station.split('.')[1]
        if loc == '':
            loc = None
        st_temp = st.select(
            station=station.split('.')[0], location=loc
        ).copy()  # [trace for trace in st if station in trace.stat.station]
        if len(st_temp) == 3:  # if len 3, put z down, rotate horizontals
            try:
                z = st_temp.select(component='Z').copy()
                st_rotated = st_rotated + z.copy()
                chans = [trace.stats.channel for trace in st_temp]
                # if not BHN BHE
                if 'H1' in str(chans) and 'H2' in str(chans):
                    try:
                        for k, trace in enumerate(st_temp):
                            if trace.stats.location == '':
                                loc = '--'
                            else:
                                loc = trace.stats.location
                            url = (
                                'http://service.iris.edu/fdsnws/station/1/query?net=%s&sta=%s&loc=%s&cha=%s&level=channel&format=text&includecomments=true&nodata=404'
                                % (
                                    trace.stats.network,
                                    trace.stats.station,
                                    loc,
                                    trace.stats.channel,
                                )
                            )
                            temp = urllib.request.urlopen(url)
                            file1 = temp.read()
                            lines = [line.split('|') for line in file1.split('\n')[1:]]
                            trace.stats.cmpaz = float(lines[0][8])
                            st_temp[k] = trace

                        z1, n1, e1 = rotate2ZNE(
                            z[0].data,
                            0,
                            -90,
                            st_temp.select(component='1')[0].data,
                            st_temp.select(component='1')[0].stats.cmpaz,
                            0,
                            st_temp.select(component='2')[0].data,
                            st_temp.select(component='2')[0].stats.cmpaz,
                            0,
                        )
                        st_temp.select(component='1')[0].data = n1
                        st_temp.select(component='1')[0].stats.channel = 'BHN'
                        st_temp.select(component='2')[0].data = e1
                        st_temp.select(component='2')[0].stats.channel = 'BHE'
                    except:
                        print(
                            'couldnt get cmpaz orientation from IRIS, rotation failed'
                        )
                        continue
                st_h = (
                    st_temp.select(component='N').copy()
                    + st_temp.select(component='E').copy()
                )
                st_h.rotate('NE->RT')
                st_rotated = st_rotated + st_h.copy()
            except:
                print('error in rotating for ' + station + ' -skipping')
        elif len(st_temp) == 1:  # if len 1, put z down, continue
            z = st_temp.select(component='Z')
            st_rotated = st_rotated + z.copy()
        elif len(st_temp) == 2:  # if len 2, probably horizontal components
            try:
                st_h = (
                    st_temp.select(component='N').copy()
                    + st_temp.select(component='E').copy()
                )
                st_h.rotate('NE->RT')
                st_rotated = st_rotated + st_h.copy()
            except:
                print(('weird number of components for ' + station + ' -skipping'))
        else:
            print(('weird number of components for ' + station + ' -skipping'))
    st_rotated.pop(0)  # pop off the placeholder

    return st_rotated


def _dip_azimuth2ZSE_base_vector(dip, azimuth):
    """
    Taken from https://github.com/obspy/obspy/blob/a8cf88bfc28b7d06b88427ca626f67418922aa52/obspy/signal/rotate.py#L188-251

    Helper function converting a vector described with azimuth and dip of unit
    length to a vector in the ZSE (Vertical, South, East) base.
    The definition of azimuth and dip is according to the SEED reference
    manual, as are the following examples (they use rounding for small
    numerical inaccuracies - also positive and negative zero are treated as
    equal):
    """

    # Convert both to radian.
    dip = np.deg2rad(dip)
    azimuth = np.deg2rad(azimuth)
    # Define the rotation axis for the dip.
    c1 = 0.0
    c2 = 0.0
    c3 = -1.0
    # Now the dip rotation matrix.
    dip_rotation_matrix = (
        np.cos(dip) * np.matrix(((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)))
        + (1 - np.cos(dip))
        * np.matrix(
            (
                (c1 * c1, c1 * c2, c1 * c3),
                (c2 * c1, c2 * c2, c2 * c3),
                (c3 * c1, c3 * c2, c3 * c3),
            )
        )
        + np.sin(dip) * np.matrix(((0, -c3, c2), (c3, 0, -c1), (-c2, c1, 0)))
    )
    # Do the same for the azimuth.
    c1 = -1.0
    c2 = 0.0
    c3 = 0.0
    azimuth_rotation_matrix = (
        np.cos(azimuth) * np.matrix(((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)))
        + (1 - np.cos(azimuth))
        * np.matrix(
            (
                (c1 * c1, c1 * c2, c1 * c3),
                (c2 * c1, c2 * c2, c2 * c3),
                (c3 * c1, c3 * c2, c3 * c3),
            )
        )
        + np.sin(azimuth) * np.matrix(((0, -c3, c2), (c3, 0, -c1), (-c2, c1, 0)))
    )
    # Now simply rotate a north pointing unit vector with both matrixes.
    temp = np.dot(azimuth_rotation_matrix, [[0.0], [-1.0], [0.0]])
    return np.array(np.dot(dip_rotation_matrix, temp)).ravel()


def readrun(filename):
    """
    Read in a saved LSforce object
    """

    with open(filename, 'rb') as f:
        result = pickle.load(f)

    return result
