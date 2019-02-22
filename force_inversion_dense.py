#import obspy.signal as obsig
import numpy as np
from obspy import Trace, Stream, read
#from scipy import sparse
import math
import glob
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import scipy as sp
import urllib.request, urllib.error, urllib.parse
import random as rnd
import pickle
from sklearn import linear_model as lm


"""
NOTE THIS HASNT YET BEEN CHECKED FOR FUNCTIONALITY IN PYTHON 3
"""


def rotate(st, baz=None):
    """
    rotate all components of st that can be rotated, to radial and transverse
    st = obspy stream object to rotate
    baz = None required if backaz already attached to st stats, list of baz corresponding to each otherwise
    """
    #implant baz in st's
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
    #get list of station location code pairs present
    staloc = list(set([trace.stats.station+'.'+trace.stats.location for trace in st]))
    st_rotated = Stream(traces=Trace())  # initialize, pop this one off later
    for station in staloc:
        #get all components with that station name
        loc = station.split('.')[1]
        if loc is '':
            loc = None
        st_temp = st.select(station=station.split('.')[0], location=loc).copy()  # [trace for trace in st if station in trace.stat.station]
        if len(st_temp) == 3:  # if len 3, put z down, rotate horizontals
            #try:
                z = st_temp.select(component='Z').copy()
                st_rotated = st_rotated+z.copy()
                chans = [trace.stats.channel for trace in st_temp]
                #if not BHN BHE
                if 'H1' in str(chans) and 'H2' in str(chans):
                    try:
                        for k, trace in enumerate(st_temp):
                            if trace.stats.location is '':
                                loc = '--'
                            else:
                                loc = trace.stats.location
                            url = ('http://service.iris.edu/fdsnws/station/1/query?net=%s&sta=%s&loc=%s&cha=%s&level=channel&format=text&includecomments=true&nodata=404' % (trace.stats.network, trace.stats.station, loc, trace.stats.channel))
                            temp = urllib.request.urlopen(url)
                            file1 = temp.read()
                            lines = [line.split('|') for line in file1.split('\n')[1:]]
                            trace.stats.cmpaz = float(lines[0][8])
                            st_temp[k] = trace

                        z1, n1, e1 = rotate2ZNE(z[0].data, 0, -90, st_temp.select(component='1')[0].data, st_temp.select(component='1')[0].stats.cmpaz, 0, st_temp.select(component='2')[0].data, st_temp.select(component='2')[0].stats.cmpaz, 0)
                        st_temp.select(component='1')[0].data = n1
                        st_temp.select(component='1')[0].stats.channel = 'BHN'
                        st_temp.select(component='2')[0].data = e1
                        st_temp.select(component='2')[0].stats.channel = 'BHE'
                    except:
                        print('couldnt get cmpaz orientation from IRIS, rotation failed')
                        continue
                st_h = st_temp.select(component='N').copy()+st_temp.select(component='E').copy()
                st_h.rotate('NE->RT')
                st_rotated = st_rotated+st_h.copy()
            #except:
                #print('error in rotating for '+station+' -skipping')
        elif len(st_temp) == 1:  # if len 1, put z down, continue
            z = st_temp.select(component='Z')
            st_rotated = st_rotated+z.copy()
        elif len(st_temp) == 2:  # if len 2, probably horizontal components
            try:
                st_h = st_temp.select(component='N').copy()+st_temp.select(component='E').copy()
                st_h.rotate('NE->RT')
                st_rotated = st_rotated+st_h.copy()
            except:
                print(('weird number of components for '+station+' -skipping'))
        else:
            print(('weird number of components for '+station+' -skipping'))
    st_rotated.pop(0)  # pop off the placeholder

    return st_rotated


def setup_timedomain(st, greendir, samplerate, weights=None, weightpre=None, period_range=[30, 100],
                     filter_order=2, zeroPhase=False, az=None):
    """
    Load in and set up matrices for time domain calculation

    Args:
        st (Stream): stream object of data to use trimmed to desired time period with horizontals
            rotated to transverse and radial and azimuths embedded
        greendir (str):
        samplerate (float): sampling rate at which to do the inversion in Hz (all data and Green's
                   functions will be resampled to this rate without filtering but after filter 
                   specified by period_range is applied)
        weights (array): array of weights to apply to each station, should be an array the same
            length as st
        weightpre (float): length of pre-noise window in seconds (if not None, noise will be used to 
                  determine weights)
        period_range (list): Range of periods to consider in inversion, in seconds
        filter_order (int): Order of filter applied over period_range
        zeroPhase (bool): If True, zeroPhase filtering will be used
        
    Returns:
        
    """
    #check if sampling rate specified is compatible with period_range
    if 1/period_range[0] > samplerate/2:
        print('samplerate and period_range are not compatible, violates Nyquist')
        return
    #filter data to band specified
    st.detrend('linear')
    st.taper(max_percentage=0.05)
    st.filter('bandpass', freqmin=1/period_range[1], freqmax=1/period_range[0], corners=filter_order, zerophase=zeroPhase)
    #resample st to samplerate
    st.resample(samplerate)
    st.taper(max_percentage=0.05)

    #make sure st data are all the same length
    lens = [len(trace.data) for trace in st]
    if len(set(lens)) != 1:
        print('traces in st are not all the same length')
        return

    K = 1e-15
    datalength = len(st[0].data)
    temp = read(glob.glob(greendir+'/*.sac')[0])[0]
    temp.resample(samplerate)
    greenlength = len(temp.data)
    if greenlength > datalength:
        print('greenlength is greater than datalength so everything is messed up')
        return
    #ADD WAY TO ACCOUNT FOR WHEN GREENLENGTH IS LONGER THAN DATALENGTH - ACTUALLY SHOULD BE AS LONG AS BOTH ADDED TOGETHER TO AVOID WRAPPING ERROR
    lenUall = datalength*len(st)
    Wvec = np.ones(lenUall)
    indx = 0
    weight = np.ones(len(st))
    for i, trace in enumerate(st):
        newline = 0
        #find component of st
        component = trace.stats.channel[2]
        station = trace.stats.station
        if component == 'Z':
            zvf = read(glob.glob(greendir+'/*'+station+'*'+'ZVF.sac')[0])[0]
            zhf = read(glob.glob(greendir+'/*'+station+'*'+'ZHF.sac')[0])[0]
            #process the same way as st
            zvf.detrend()
            zvf.taper(max_percentage=0.05)
            zvf.filter('bandpass', freqmin=1/period_range[1], freqmax=1/period_range[0],
                       corners=filter_order, zerophase=zeroPhase)
            zvf.resample(samplerate)
            zvf.taper(max_percentage=0.05)

            zhf.detrend()
            zhf.taper(max_percentage=0.05)
            zhf.filter('bandpass', freqmin=1/period_range[1], freqmax=1/period_range[0],
                       corners=filter_order, zerophase=zeroPhase)
            zhf.resample(samplerate)
            zhf.taper(max_percentage=0.05)

            ZVF = makeconvmat(zvf.data, size=(datalength, datalength))  # sparse.diags(zvff,0)
            ZHF = makeconvmat(zhf.data, size=(datalength, datalength))  # sparse.diags(zhff,0)
            az = math.radians(trace.stats.azimuth)
            newline = np.hstack((K*ZVF, K*ZHF*math.cos(az), K*ZHF*math.sin(az)))  # sparse.hstack((K*ZVF, K*ZHF*math.cos(az), K*ZHF*math.sin(az)))
            datline = trace.data
        elif component == 'R':
            #import pdb; pdb.set_trace()
            rvf = read(glob.glob(greendir+'/*'+station+'*'+'RVF.sac')[0])[0]
            rhf = read(glob.glob(greendir+'/*'+station+'*'+'RHF.sac')[0])[0]
            #process the same way as st
            rvf.detrend()
            rvf.taper(max_percentage=0.05)

            rvf.filter('bandpass', freqmin=1/period_range[1], freqmax=1/period_range[0],
                       corners=filter_order, zerophase=zeroPhase)
            rvf.resample(samplerate)
            rvf.taper(max_percentage=0.05)

            rhf.detrend()
            rhf.taper(max_percentage=0.05)
            rhf.filter('bandpass', freqmin=1/period_range[1], freqmax=1/period_range[0],
                       corners=filter_order, zerophase=zeroPhase)
            rhf.resample(samplerate)
            rhf.taper(max_percentage=0.05)

            RVF = makeconvmat(rvf.data, size=(datalength, datalength))  # sparse.diags(rvff,0)
            RHF = makeconvmat(rhf.data, size=(datalength, datalength))  # sparse.diags(rhff,0)
            az = math.radians(trace.stats.azimuth)
            newline = np.hstack((K*RVF, K*RHF*math.cos(az), K*RHF*math.sin(az)))  # sparse.hstack((K*RVF, K*RHF*math.cos(az), K*RHF*math.sin(az)))
            datline = trace.data
        elif component == 'T':
            thf = read(glob.glob(greendir+'/*'+station+'*'+'THF.sac')[0])[0]
            #process the same way as st
            thf.detrend()
            thf.taper(max_percentage=0.05)
            thf.filter('bandpass', freqmin=1/period_range[1], freqmax=1/period_range[0],
                       corners=filter_order, zerophase=zeroPhase)
            thf.resample(samplerate)
            thf.taper(max_percentage=0.05)
            THF = makeconvmat(thf.data, size=(datalength, datalength))  # sparse.diags(thff,0)
            TVF = 0*THF.copy()
            az = math.radians(trace.stats.azimuth)
            newline = np.hstack((TVF, K*THF*math.sin(az), -K*THF*math.cos(az)))  # sparse.hstack((K*TVF, K*THF*math.sin(az), -K*THF*math.cos(az)))
            datline = trace.data
        else:
            print('st not rotated to T and R, abort!')
            import pdb
            pdb.set_trace()
            return
        if i == 0:
            G = newline.copy()
            d = datline.copy()
        else:
            G = np.vstack((G, newline.copy()))  # sparse.vstack((G,newline))
            d = np.hstack((d, datline.copy()))
        if weights is not None:
            if weights is 'prenoise':  #CHANGE THIS TO USE STD?
                weight[i] = 10**-7*(1./np.mean(np.abs(trace.data[0:int(weightpre*trace.stats.sampling_rate)])))
            elif weights is 'distance':
                weight[i] = trace.stats.rdist
            elif len(weights) > 1:
                weight[i] = weights[i]
            Wvec[indx:indx+datalength] = Wvec[indx:indx+datalength]*weight[i]
            indx += datalength

    if np.shape(G)[0] != len(d):
        print('G and d sizes are not compatible, fix something somewhere')
    G = G * 1/samplerate  # multiply by sample interval (sec) since convolution is an integral
    d = d * 100  # convert data from m to cm
    if weights is not None:
        W = np.diag(Wvec)  # sparse.diags(Wvec,0)
    else:
        W = None
    #import pdb; pdb.set_trace()
    return G, d, W, weight


def setup_freqdomain(st, greendir, samplerate, weights=None, weightpre=None, period_range=[30, 100],
                     filter_order=2, zeroPhase=False, az=None):
    """
    Load in and set up matrices for long-period landslide inversion in the frequency domain
    
    Args:
        st (Stream): stream object of data to use trimmed to desired time period with horizontals
            rotated to transverse and radial and azimuths embedded
        greendir (str):
        samplerate (float): sampling rate at which to do the inversion in Hz (all data and Green's
                   functions will be resampled to this rate without filtering but after filter 
                   specified by period_range is applied)
        weights (array): array of weights to apply to each station, should be an array the same
            length as st
        weightpre (float): length of pre-noise window in seconds (if not None, noise will be used to 
                  determine weights)
        period_range (list): Range of periods to consider in inversion, in seconds
        filter_order (int): Order of filter applied over period_range
        zeroPhase (bool): If True, zeroPhase filtering will be used
        
    Returns:
        
    """
    #check if sampling rate specified is compatible with period_range
    if 1/period_range[0] > samplerate/2:
        print('samplerate and period_range are not compatible, violates Nyquist')
        return
    #filter data to band specified
    st.detrend('linear')
    st.taper(max_percentage=0.05)
    st.filter('bandpass', freqmin=1/period_range[1], freqmax=1/period_range[0],
              corners=filter_order, zerophase=zeroPhase)
    # resample st to samplerate
    st.resample(samplerate)
    st.taper(max_percentage=0.05)

    # make sure st data are all the same length
    lens = [len(trace.data) for trace in st]
    if len(set(lens)) != 1:
        print('traces in st are not all the same length')
        return

    K = 1e-15
    datalength = len(st[0].data)
    temp = read(glob.glob(greendir+'/*.sac')[0])[0]
    temp.resample(samplerate)
    #greenlength = len(temp.data)
    NFFT = nextpow2(datalength)  # +greenlength) #needs to be the length of the two added together because convolution length M+N-1
    lenUall = NFFT*len(st)
    #if NFFT%2==0: #this is only for rfft
        #lenUall = (NFFT/2+1)*len(st)
    #else:
        #lenUall = (NFFT+1)/2*len(st)
    Wvec = np.ones(lenUall)
    indx = 0
    weight = np.ones(len(st))
    for i, trace in enumerate(st):
        #find component of st
        component = trace.stats.channel[2]
        station = trace.stats.station
        if component == 'Z':
            zvf = read(glob.glob(greendir+'/*'+station+'*'+'ZVF.sac')[0])[0]
            zhf = read(glob.glob(greendir+'/*'+station+'*'+'ZHF.sac')[0])[0]
            #process the same way as st
            zvf.taper(max_percentage=0.05)
            zvf.filter('bandpass', freqmin=1/period_range[1], freqmax=1/period_range[0],
                       corners=filter_order, zerophase=zeroPhase)
            zvf.resample(samplerate)
            zhf.taper(max_percentage=0.05)
            zhf.filter('bandpass', freqmin=1/period_range[1], freqmax=1/period_range[0],
                       corners=filter_order, zerophase=zeroPhase)
            zhf.resample(samplerate)
            zvff = np.fft.fft(zvf.data, NFFT)
            zhff = np.fft.fft(zhf.data, NFFT)
            ZVF = np.diag(zvff)  # sparse.diags(zvff,0)
            ZHF = np.diag(zhff)  # sparse.diags(zhff,0)
            az = math.radians(np.round(trace.stats.azimuth))
            newline = np.hstack((K*ZVF, K*ZHF*math.cos(az), K*ZHF*math.sin(az)))  # sparse.hstack((K*ZVF, K*ZHF*math.cos(az), K*ZHF*math.sin(az)))
        elif component == 'R':
            rvf = read(glob.glob(greendir+'/*'+station+'*'+'RVF.sac')[0])[0]
            rhf = read(glob.glob(greendir+'/*'+station+'*'+'RHF.sac')[0])[0]
            #process the same way as st
            rvf.taper(max_percentage=0.05)
            rvf.filter('bandpass', freqmin=1/period_range[1], freqmax=1/period_range[0], corners=filter_order, zerophase=zeroPhase)
            rvf.resample(samplerate)
            rhf.taper(max_percentage=0.05)
            rhf.filter('bandpass', freqmin=1/period_range[1], freqmax=1/period_range[0],
                       corners=filter_order, zerophase=zeroPhase)
            rhf.resample(samplerate)
            rvff = np.fft.fft(rvf.data, NFFT)
            rhff = np.fft.fft(rhf.data, NFFT)
            RVF = np.diag(rvff)  # sparse.diags(rvff,0)
            RHF = np.diag(rhff)  # sparse.diags(rhff,0)
            az = math.radians(np.round(trace.stats.azimuth))
            newline = np.hstack((K*RVF, K*RHF*math.cos(az), K*RHF*math.sin(az)))  # sparse.hstack((K*RVF, K*RHF*math.cos(az), K*RHF*math.sin(az)))
        elif component == 'T':
            thf = read(glob.glob(greendir+'/*'+station+'*'+'THF.sac')[0])[0]
            #process the same way as st
            thf.taper(max_percentage=0.05)
            thf.filter('bandpass', freqmin=1/period_range[1], freqmax=1/period_range[0],
                       corners=filter_order, zerophase=zeroPhase)
            thf.resample(samplerate)
            thff = np.fft.fft(thf.data, NFFT)
            THF = np.diag(thff)  # sparse.diags(thff,0)
            TVF = 0*THF
            az = math.radians(np.round(trace.stats.azimuth))
            newline = np.hstack((TVF, K*THF*math.sin(az), -K*THF*math.cos(az)))  # sparse.hstack((K*TVF, K*THF*math.sin(az), -K*THF*math.cos(az)))
        else:
            print('st not rotated to T and R, abort!')
            return
        datline = np.fft.fft(trace.data, NFFT)
        if i == 0:
            G = newline.copy()
            d = datline.copy()
        else:
            G = np.vstack((G, newline.copy()))  # sparse.vstack((G,newline))
            d = np.hstack((d, datline.copy()))
        if weights is not None:
            if weights is 'prenoise':
                weight[i] = 10**-7*(1./np.mean(np.abs(trace.data[0:int(weightpre*trace.stats.sampling_rate)])))
            elif weights is 'distance':
                weight[i] = trace.stats.rdist
            elif len(weights) > 1:
                weight[i] = weights[i]
            Wvec[indx:indx+NFFT] = Wvec[indx:indx+NFFT]*weight[i]
            indx += NFFT
    if np.shape(G)[0] != len(d):
        print('G and d sizes are not compatible, fix something somewhere')
        return
    G = G * 1./samplerate  # multiply by sample interval (sec) since convolution is an integral
    d = d * 100  # convert data from m to cm
    if weights is not None:
        W = np.diag(Wvec)  # sparse.diags(Wvec,0)
    else:
        W = None

    return G, d, W


def invert(G, d, samplerate, numsta, datlenorig, W=None, T0=0, Tikhratio=[1.0, 0., 0.], domain='time',
           alphaset=None, zeroTime=None, imposeZero=False, addtoZero=False, alpha_method='Lcurve',
           maxduration=None):
    """
    Wrapper function to perform single force inversion of long-period landslide seismic signal

    Args:
        G (array): model matrix (m x n)
        d (array): vector of concatenated data (m x 1)
        samplerate (float): samplerate in Hz of seismic data and green's functions (the two must be equal)
        numsta (int): number of channels used
        datlenorig (int): length, in samples, or original data prior to zero padding etc.
        W: optional weighting matrix, same size as G
        T0 (float): Reference time T0 used in Green's function computation (usually 0.)
        Tikhratio (array): Proportion each regularization method contributes, where values correspond
            to [zeroth, first order, second order]. Must add to 1. 
        domain (str): specifies whether calculations are in time domain ('time', default) or
            freq domain ('freq')
        alphaset (float): Set regularization parameter, if None, will search for best alpha
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

    if np.sum(Tikhratio) != 1.:
        raise Exception('Tikhonov ratios must add to 1')

    if W is not None:
        Ghat = W.dot(G)  # np.dot(W.tocsr(),G.tocsr())
        dhat = W.dot(d)  # np.dot(W.tocsr(),sparse.csr_matrix(d))
    else:
        Ghat = G  # G.tocsr()
        dhat = d  # sparse.csr_matrix(d)

    m, n = np.shape(Ghat)

    Ghatnorm = np.linalg.norm(Ghat)
    #Ghatmax = np.abs(Ghat).max()
    
    if addtoZero is True:  # constrain forces to add to zero
        scaler = Ghatnorm #10**(np.round(np.log10(Ghatmax)+4))#10**(np.round(np.log10(Ghatnorm))) # 10**(np.round(np.log10(Ghatmax)+4))
        first1 = np.hstack((np.ones(datlenorig), np.zeros(2*datlenorig)))
        second1 = np.hstack((np.zeros(datlenorig), np.ones(datlenorig), np.zeros(datlenorig)))
        third1 = np.hstack((np.zeros(2*datlenorig), np.ones(datlenorig)))
        A = np.vstack((first1, second1, third1))*scaler
        Ghat = np.vstack((Ghat, A))
        dhat = np.hstack((dhat, np.zeros(3)))
        #import pdb; pdb.set_trace()

    scaler = Ghatnorm/15.#10**(np.round(np.log10(Ghatmax))+0.5) #10**(np.round(np.log10((Ghatnorm)+0.5)))
    if imposeZero is True:  # tell model when there should be no forces
        if zeroTime is None:
            raise Exception('imposeZero set to True but no zeroTime provided')
        len2 = int(np.round((zeroTime)*samplerate))
        len3 = int(np.round(0.2*len2))  # 20% taper overlapping into main event by x seconds
        #halflen3 = int(len3/2)
        temp = np.hanning(2*len3)
        temp = temp[len3:]
        vals = np.hstack((np.ones(len2-len3), temp))
        for i, val in enumerate(vals):
            first1 = np.zeros(3*datlenorig)
            second1 = first1.copy()
            third1 = first1.copy()
            first1[i] = val
            second1[i+datlenorig] = val
            third1[i+2*datlenorig] = val
            if i == 0:
                A = np.vstack((first1, second1, third1))
            else:
                A = np.vstack((A, first1, second1, third1))
        A = A*scaler
        Ghat = np.vstack((Ghat, A))
        dhat = np.hstack((dhat, np.zeros(len(vals)*3)))
        
    if maxduration is not None:
        if zeroTime is None:
            zeroTime = 0.
        startind = int((zeroTime + maxduration)*samplerate)
        len2 = int(datlenorig - startind)
        len3 = int(np.round(0.2*len2))  # 20% taper so zero imposition isn't sudden
        temp = np.hanning(2*len3)
        temp = temp[:len3]
        vals = np.hstack((temp, np.ones(len2-len3)))
        for i, val in enumerate(vals):
            place = i + startind
            first1 = np.zeros(3*datlenorig)
            second1 = first1.copy()
            third1 = first1.copy()
            first1[place] = val
            second1[place+datlenorig] = val
            third1[place+2*datlenorig] = val
            if i == 0:
                A = np.vstack((first1, second1, third1))
            else:
                A = np.vstack((A, first1, second1, third1))
        A = A*scaler
        Ghat = np.vstack((Ghat, A))
        dhat = np.hstack((dhat, np.zeros(len(vals)*3)))

    if alphaset is not None:
        alpha = alphaset
    dhat = dhat.T
    
    # Build roughening matrix
    I = np.eye(n, n) # sparse.eye(np.shape(G)[1],np.shape(G)[1])
    if Tikhratio[1] != 0.:
        # Build L1 (first order) roughening matrix
        L1 = np.diag(-1 * np.ones(n)) + np.diag(np.ones(n-1), k=1)
        L1part = np.dot(L1.T, L1)
    else:
        L1part = 0.
        L1 = 0.
    if Tikhratio[2] != 0.:
        # Build L2 (second order) roughening matrix
        L2 = np.diag(np.ones(n)) + np.diag(-2*np.ones(n-1), k=1) + np.diag(np.ones(n-2), k=2)
        L2part = np.dot(L2.T, L2)
    else:
        L2 = 0.
        L2part = 0.

    if alphaset is None:
        if alpha_method == 'Lcurve':
            alpha, fit1, size1, alphas = findalpha(Ghat, dhat, I, L1, L2, Tikhratio)
        else:
            alpha, fit1, size1, alphas = findalphaD(Ghat, dhat, I, zeroTime,
                                                    samplerate, numsta, datlenorig, L1=L1, L2=L2,
                                                    Tikhratio=Tikhratio)#, tolerance=0.5)
        print('best alpha is %6.1e' % alpha)
    else:
        fit1 = None
        size1 = None
        alphas = None

    Ghat = np.matrix(Ghat)
    Apart = np.dot(Ghat.H, Ghat)
    
    A = Apart+alpha**2*(Tikhratio[0]*I + Tikhratio[1]*L1part + Tikhratio[2]*L2part) # Combo of all regularization things (if any are zero they won't matter)    
    x = np.squeeze(np.asarray(np.dot(Ghat.H, dhat)))

    if domain is 'freq':
        model, residuals, rank, s = sp.linalg.lstsq(A, x)  # sparse.linalg.spsolve(Ghat.T*Ghat+alpha**2*I,Ghat.T*dhat)
        #import pdb;pdb.set_trace()
        div = len(model)/3
        Zforce = - np.real(np.fft.ifft(model[0:div])/10**5)  # convert from dynes to newtons, flip so up is positive
        Nforce = np.real(np.fft.ifft(model[div:2*div])/10**5)
        Eforce = np.real(np.fft.ifft(model[2*div:])/10**5)
        #run forward model
        df_new = np.dot(G, model.T)  # forward_model(G,model)
        #convert d and df_new back to time domain
        dt, dtnew = back2time(d, df_new, numsta, datlenorig)

    else:  # domain is time
        model, residuals, rank, s = sp.linalg.lstsq(A, x)
        #model,residuals,rank,s = sp.linalg.lstsq(Ghat,dhat,cond=alpha)
        div = int(len(model)/3)
        Zforce = - model[0:div]/10**5  # convert from dynes to netwons, flip so up is positive
        Nforce = model[div:2*div]/10**5
        Eforce = model[2*div:]/10**5
        dtnew = G.dot(model)  # forward_model(G,model)
        dtnew = np.reshape(dtnew, (numsta, datlenorig))
        dt = np.reshape(d, (numsta, datlenorig))

    #compute variance reduction
    VR = varred(dt, dtnew)
    print(('variance reduction %f percent') % (VR,))
    tvec = np.arange(0, len(Zforce)*1/samplerate, 1/samplerate)-T0
    #np.linspace(0,(len(Zforce)-1)*1/samplerate,len(Zforce))-T0
    if zeroTime is not None:
        tvec = tvec - zeroTime
    return model, Zforce, Nforce, Eforce, tvec, VR, dt, dtnew, alpha, fit1, size1, alphas


def invert_lasso(G, d, samplerate, numsta, datlenorig, W=None, T0=0, alpharatio=10., domain='time',
                 alphaset=None, zeroTime=None, imposeZero=False, addtoZero=False, alpha_method='Lcurve',
                 maxduration=None):
    """
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

    if W is not None:
        Ghat = W.dot(G)  # np.dot(W.tocsr(),G.tocsr())
        dhat = W.dot(d)  # np.dot(W.tocsr(),sparse.csr_matrix(d))
    else:
        Ghat = G  # G.tocsr()
        dhat = d  # sparse.csr_matrix(d)

    m, n = np.shape(Ghat)

    Ghatnorm = np.linalg.norm(Ghat)
    #Ghatmax = np.abs(Ghat).max()
    
    if addtoZero is True:  # constrain forces to add to zero
        scaler = Ghatnorm #10**(np.round(np.log10(Ghatmax)+4))#10**(np.round(np.log10(Ghatnorm))) # 10**(np.round(np.log10(Ghatmax)+4))
        first1 = np.hstack((np.ones(datlenorig), np.zeros(2*datlenorig)))
        second1 = np.hstack((np.zeros(datlenorig), np.ones(datlenorig), np.zeros(datlenorig)))
        third1 = np.hstack((np.zeros(2*datlenorig), np.ones(datlenorig)))
        A = np.vstack((first1, second1, third1))*scaler
        Ghat = np.vstack((Ghat, A))
        dhat = np.hstack((dhat, np.zeros(3)))
        #import pdb; pdb.set_trace()

    scaler = Ghatnorm/15.#10**(np.round(np.log10(Ghatmax))+0.5) #10**(np.round(np.log10((Ghatnorm)+0.5)))
    if imposeZero is True:  # tell model when there should be no forces
        if zeroTime is None:
            raise Exception('imposeZero set to True but no zeroTime provided')
        len2 = int(np.round((zeroTime)*samplerate))
        len3 = int(np.round(0.2*len2))  # 20% taper overlapping into main event by x seconds
        #halflen3 = int(len3/2)
        temp = np.hanning(2*len3)
        temp = temp[len3:]
        vals = np.hstack((np.ones(len2-len3), temp))
        for i, val in enumerate(vals):
            first1 = np.zeros(3*datlenorig)
            second1 = first1.copy()
            third1 = first1.copy()
            first1[i] = val
            second1[i+datlenorig] = val
            third1[i+2*datlenorig] = val
            if i == 0:
                A = np.vstack((first1, second1, third1))
            else:
                A = np.vstack((A, first1, second1, third1))
        A = A*scaler
        Ghat = np.vstack((Ghat, A))
        dhat = np.hstack((dhat, np.zeros(len(vals)*3)))
        
    if maxduration is not None:
        if zeroTime is None:
            zeroTime = 0.
        startind = int((zeroTime + maxduration)*samplerate)
        len2 = int(datlenorig - startind)
        len3 = int(np.round(0.2*len2))  # 20% taper so zero imposition isn't sudden
        temp = np.hanning(2*len3)
        temp = temp[:len3]
        vals = np.hstack((temp, np.ones(len2-len3)))
        for i, val in enumerate(vals):
            place = i + startind
            first1 = np.zeros(3*datlenorig)
            second1 = first1.copy()
            third1 = first1.copy()
            first1[place] = val
            second1[place+datlenorig] = val
            third1[place+2*datlenorig] = val
            if i == 0:
                A = np.vstack((first1, second1, third1))
            else:
                A = np.vstack((A, first1, second1, third1))
        A = A*scaler
        Ghat = np.vstack((Ghat, A))
        dhat = np.hstack((dhat, np.zeros(len(vals)*3)))

    if alphaset is not None:
        alpha = alphaset
    dhat = dhat.T
    
    # Build roughening matrix
    if alpharatio is not None:
        # Build L2 (second order) roughening matrix
        L2 = np.diag(np.ones(n)) + np.diag(-2*np.ones(n-1), k=1) + np.diag(np.ones(n-2), k=2)

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
            raise Exception('alphaset=None not implemented yet for smoothed Lasso. Must assign alpha')
            #print('best alpha is %6.1e' % alpha)
    else:
        if alpharatio is not None:
            Ghat2 = np.vstack((Ghat, alphaset*L2))
            dhat2 = np.hstack((dhat, np.zeros(np.shape(L2[0]))))
        else:
            Ghat2 = Ghat
            dhat2 = dhat
            alpharatio = 1.
        lasso = lm.Lasso(alpha=(alphaset*alpharatio)**2)
        lasso.fit(Ghat2, dhat2)
        fit1 = None
        size1 = None
        alphas = None

    model = lasso.coef_
    div = int(len(model)/3)

    if domain is 'freq':
        Zforce = - np.real(np.fft.ifft(model[0:div])/10**5)  # convert from dynes to newtons, flip so up is positive
        Nforce = np.real(np.fft.ifft(model[div:2*div])/10**5)
        Eforce = np.real(np.fft.ifft(model[2*div:])/10**5)
        #run forward model
        df_new = np.dot(G, model.T)  # forward_model(G,model)
        #convert d and df_new back to time domain
        dt, dtnew = back2time(d, df_new, numsta, datlenorig)

    else:  # domain is time
        Zforce = - model[0:div]/10**5  # convert from dynes to netwons, flip so up is positive
        Nforce = model[div:2*div]/10**5
        Eforce = model[2*div:]/10**5
        dtnew = G.dot(model)  # forward_model(G,model)
        dtnew = np.reshape(dtnew, (numsta, datlenorig))
        dt = np.reshape(d, (numsta, datlenorig))

    #compute variance reduction
    VR = varred(dt, dtnew)
    print(('variance reduction %f percent') % (VR,))
    tvec = np.arange(0, len(Zforce)*1/samplerate, 1/samplerate)-T0
    #np.linspace(0,(len(Zforce)-1)*1/samplerate,len(Zforce))-T0
    if zeroTime is not None:
        tvec = tvec - zeroTime
    return model, Zforce, Nforce, Eforce, tvec, VR, dt, dtnew, alpha, fit1, size1, alphas


def jackknife(G, d, samplerate, numsta, datlenorig, num_iter=200, frac_delete=0.5, W=None, T0=0,
              L0L2ratio=[0.9, 0.1], domain='time', alphaset=None, zeroTime=None, imposeZero=False,
              addtoZero=False):
    # NEED TO UPDATE FOR CONSISTENCY WITH invert
    #inversion with jackknife to estimate uncertainties due to data selection
    #NEED TO SET ALPHA
    #First iteration is full inversion - POP IT OUT AT THE END
    #domain specifies whether calculations are in time domain (default) or freq domain
    #L1 and L2 norms, include as % of each - NEED TO ADD THIS STILL
    #addtozero = make forces add to zero

    # Initialize some things
    datlen = len(d)/numsta
    Zforce_all = []
    ZforceL = []
    ZforceU = []
    Nforce_all = []
    NforceL = []
    NforceU = []
    Eforce_all = []
    EforceL = []
    EforceU = []
    VR_all = []

    if alphaset is None:
        print('You need to define alpha')
    else:
        alpha = alphaset

    # Apply weighting to matrices
    if W is not None:
        Ghat = W.dot(G)  # np.dot(W.tocsr(),G.tocsr())
        dhat = W.dot(d)  # np.dot(W.tocsr(),sparse.csr_matrix(d))
    else:
        Ghat = G  # G.tocsr()
        dhat = d  # sparse.csr_matrix(d)

    Ghatmeanfull = np.abs(Ghat).mean()

    # Start jackknife iterations
    for ii in range(num_iter):
        if ii == 0:  # first time do full inversion
            Ghat1 = Ghat
            dhat1 = dhat
        else:
            numcut = int(round(frac_delete*numsta))
            indxcut = rnd.sample(list(range(numsta)), numcut)
            #TEST AND MAKE SURE THIS IS RIGHT
            obj = [sum(ind) for ind in zip(np.tile(list(range(datlen)), len(indxcut)), np.repeat([x*datlen for x in indxcut], datlen))]
            dhat1 = np.delete(dhat.copy(), obj)
            Ghat1 = np.delete(Ghat.copy(), obj, axis=0)

        Ghatmax = np.abs(Ghat1).max()
        Ghatmean = np.abs(Ghat1).mean()
        #print ('Ghatmax is %s' % (Ghatmax,))
        #print ('Ghatmean is %s' % (Ghatmean,))

        if addtoZero is True:  # constrain forces to add to zero
            scaler = 10**(np.round(np.log10(Ghatmax)+4))
            first1 = np.hstack((np.ones(datlenorig), np.zeros(2*datlenorig)))
            second1 = np.hstack((np.zeros(datlenorig), np.ones(datlenorig), np.zeros(datlenorig)))
            third1 = np.hstack((np.zeros(2*datlenorig), np.ones(datlenorig)))
            A = np.vstack((first1, second1, third1))*scaler
            Ghat1 = np.vstack((Ghat1, A))
            dhat1 = np.hstack((dhat1, np.zeros(3)))

        if imposeZero is True:  # tell model when there should be no forces
            scaler = 10**(np.round(np.log10(Ghatmax))+0.5)
            #import pdb; pdb.set_trace()
            len2 = int(np.round((zeroTime)*samplerate))
            len3 = int(np.round(0.2*len2))  # 20% taper overlapping into main event by x(0) seconds
            temp = np.hanning(2*len3)
            temp = temp[len3:]
            vals = np.hstack((np.ones(len2-len3), temp))
            for i, val in enumerate(vals):
                first1 = np.zeros(3*datlenorig)
                second1 = first1.copy()
                third1 = first1.copy()
                first1[i] = val
                second1[i+datlenorig] = val
                third1[i+2*datlenorig] = val
                if i == 0:
                    A = np.vstack((first1, second1, third1))
                else:
                    A = np.vstack((A, first1, second1, third1))
            A = A*scaler
            Ghat1 = np.vstack((Ghat1, A))
            dhat1 = np.hstack((dhat1, np.zeros(len(vals)*3)))

        dhat = dhat.T
        I = np.eye(np.shape(Ghat)[1], np.shape(Ghat)[1])  # sparse.eye(np.shape(G)[1],np.shape(G)[1])

        alphatemp = Ghatmean/Ghatmeanfull*alpha
        # print('alphatemp is %s' % (alphatemp,))

        if domain is 'freq':
            model, residuals, rank, s = sp.linalg.lstsq(np.dot(Ghat1.T, Ghat1)+alphatemp**2*I, np.dot(Ghat1.T, dhat1))  # sparse.linalg.spsolve(Ghat.T*Ghat+alpha**2*I,Ghat.T*dhat)
            div = len(model)/3
            Zforce = - np.real(np.fft.ifft(model[0:div])/10**5)  # convert from dynes to newtons, flip so up is positive
            Nforce = np.real(np.fft.ifft(model[div:2*div])/10**5)
            Eforce = np.real(np.fft.ifft(model[2*div:])/10**5)
            #run forward model
            df_new = np.dot(G, model.T)  # forward_model(G,model)
            #convert d and df_new back to time domain
            dt, dtnew = back2time(d, df_new, numsta, datlenorig)

        else:  # domain is time
            model, residuals, rank, s = sp.linalg.lstsq(np.dot(Ghat1.T, Ghat1)+alphatemp**2*I, np.dot(Ghat1.T, dhat1))
            div = len(model)/3
            Zforce = - model[0:div]/10**5  # convert from dynes to netwons, flip so up is positive
            Nforce = model[div:2*div]/10**5
            Eforce = model[2*div:]/10**5
            dtnew = G.dot(model)  # forward_model(G,model)
            dtnew = np.reshape(dtnew, (numsta, datlenorig))
            dt = np.reshape(d, (numsta, datlenorig))

        #compute variance reduction
        VR = varred(dt, dtnew)
        print(('variance reduction %f percent') % (VR,))

        # Save all critical information
        if ii == 0:  # distinguish data from full run
            model_full = model.copy()
            Zforce_full = Zforce.copy()
            Nforce_full = Nforce.copy()
            Eforce_full = Eforce.copy()
            VR_full = VR.copy()
            dtnew_full = dtnew.copy()
        else:  # Save run info in structure of some sort
            Zforce_all.append(Zforce.copy())
            Nforce_all.append(Nforce.copy())
            Eforce_all.append(Eforce.copy())
            VR_all.append(VR.copy())
            #dtnew_all.append(dtnew.copy())
            #alpha_all.append(alpha.copy())

    # Pull out upper and lower limits (middle 95%) and save as two lists
    for ll in range(len(Zforce_all[0])):
        tempZ = [val[ll] for val in Zforce_all]
        tempZ.sort()
        ZforceL.append(tempZ[int(round(0.025*len(tempZ)))]-1)  # subtract one because indexing starts at 0
        ZforceU.append(tempZ[int(round(0.975*len(tempZ)))-1])
        tempE = [val[ll] for val in Eforce_all]
        tempE.sort()
        EforceL.append(tempE[int(round(0.025*len(tempE)))]-1)
        EforceU.append(tempE[int(round(0.975*len(tempE)))-1])
        tempN = [val[ll] for val in Nforce_all]
        tempN.sort()
        NforceL.append(tempN[int(round(0.025*len(tempN)))]-1)
        NforceU.append(tempN[int(round(0.975*len(tempN)))-1])

    ZforceL = np.array(ZforceL)
    ZforceU = np.array(ZforceU)
    NforceL = np.array(NforceL)
    NforceU = np.array(NforceU)
    EforceL = np.array(EforceL)
    EforceU = np.array(EforceU)

    tvec = np.arange(0, len(Zforce_full)*1/samplerate, 1/samplerate)-T0
    #np.linspace(0,(len(Zforce)-1)*1/samplerate,len(Zforce))-T0
    if zeroTime is not None:
        tvec = tvec - zeroTime

    return model_full, Zforce_full, ZforceL, ZforceU, Nforce_full, NforceL, NforceU, Eforce_full, EforceL, EforceU, tvec, VR_full, dt, dtnew_full  # , alpha , fit1, size1, alphas, curves


def findalpha(Ghat, dhat, I, L1=0., L2=0., Tikhratio=[1., 0., 0.]):
    """
    Find best regularization (trade-off) parameter, alpha, by computing model with many values of
    alpha, plotting Lcurve, and finding point of steepest curvature where slope is negative.
    
    Args:
        Ghat (array): m x n matrix of 
        dhat (array): 1 x n array of weighted data
        I (array): Identity matrix
        L1 (array): First order roughening matrix, if 0., will use only zeroth order Tikhonov reg.
        L2 (array): Second order roughening matrix, if 0., will use only zeroth order Tikhonov reg.
        Tikhratio (list): Proportion each regularization method contributes, where values correspond
            to [zeroth, first order, second order]. Must add to 1.
    Returns:
        
    """
    templ1 = np.ceil(np.log10(np.linalg.norm(Ghat)))  # Roughly estimate largest singular value (should not use alpha larger than expected largest singular value)
    templ2 = np.arange(templ1-6, templ1-2)
    alphas = 10**templ2
    fit1 = []
    size1 = []
    
    # Convert to matrix so can use complex conjugate .H (so this code will work for both time and freq domains)
    Ghat = np.matrix(Ghat)

    Apart = np.dot(Ghat.H, Ghat)
    if type(L2) == float or type(L2) == int:
        L2part = 0.
    else:
        L2part = np.dot(L2.T, L2)
    if type(L1) == float or type(L1) == int:
        L1part = 0.
    else:
        L1part = np.dot(L1.T, L1)

    x = np.squeeze(np.asarray(np.dot(Ghat.H, dhat)))
    
    #rough first iteration
    for alpha in alphas:
        A = Apart+alpha**2*(Tikhratio[0]*I + Tikhratio[1]*L1part + Tikhratio[2]*L2part) # Combo of all regularization things
        model, residuals, rank, s = sp.linalg.lstsq(A, x)  # sparse.csr_matrix(sparse.linalg.spsolve(Ghat.T*Ghat+alpha**2*I,Ghat.T*dhat))
        temp1 = np.dot(Ghat, model.T)-dhat  # np.dot(Ghat.todense(),model.todense().T)-dhat.todense()
        fit1.append(sp.linalg.norm(temp1))
        size1.append(sp.linalg.norm(Tikhratio[0]*model) + sp.linalg.norm(Tikhratio[1]*np.dot(L1part, model))\
                     + sp.linalg.norm(Tikhratio[2]*np.dot(L2part, model)))  # size1.append(sp.linalg.norm(model.todense()))
    fit1 = np.array(fit1)
    size1 = np.array(size1)
    curves = curvature(np.log10(fit1), np.log10(size1))
    # Zero out any points where function is concave so avoid picking points form dropoff at end
    slp2 = np.gradient(np.gradient(np.log10(size1), np.log10(fit1)), np.log10(fit1))
    alphas = np.array(alphas)
    tempcurve = curves.copy()
    tempcurve[slp2 < 0] = np.max(curves)
    idx = np.argmin(tempcurve)
    alpha = alphas[idx] #[alpha for i, alpha in enumerate(alphas) if curves[i] == curves.min()]

    # Then hone in
    alphas = np.logspace(np.round(np.log10(alpha))-1, np.round(np.log10(alpha))+1, 20)
    fit1 = []
    size1 = []
    for newalpha in alphas:
        A = Apart+newalpha**2*(Tikhratio[0]*I + Tikhratio[1]*L1part + Tikhratio[2]*L2part) # Combo of all regularization things
        model, residuals, rank, s = sp.linalg.lstsq(A, x) # sparse.csr_matrix(sparse.linalg.spsolve(Ghat.T*Ghat+alpha**2*I,Ghat.T*dhat))
        temp1 = np.dot(Ghat, model.T)-dhat  # np.dot(Ghat.todense(),model.todense().T)-dhat.todense()
        fit1.append(sp.linalg.norm(temp1))
        size1.append(sp.linalg.norm(Tikhratio[0]*model) + sp.linalg.norm(Tikhratio[1]*np.dot(L1part, model))\
                     + sp.linalg.norm(Tikhratio[2]*np.dot(L2part, model)))  # size1.append(sp.linalg.norm(model.todense()))
        #size1.append(sp.linalg.norm(Tikhratio[0]*model + Tikhratio[1]*L1part*model + Tikhratio[2]*L2part*model))  # size1.append(sp.linalg.norm(model.todense()))
    fit1 = np.array(fit1)
    size1 = np.array(size1)
    curves = curvature(np.log10(fit1), np.log10(size1))
    # Zero out any points where function is concave so avoid picking points form dropoff at end
    slp2 = np.gradient(np.gradient(np.log10(size1), np.log10(fit1)), np.log10(fit1))
    alphas = np.array(alphas)
    tempcurve = curves.copy()
    tempcurve[slp2 < 0] = np.max(curves)
    #import pdb; pdb.set_trace()
    #curves[curves < 1]
    idx = np.argmin(tempcurve)
    bestalpha = alphas[idx]#[alpha1 for i, alpha1 in enumerate(alphas) if tempcurves[i] == curves.min()]

    Lcurve(fit1, size1, alphas)
    if type(bestalpha) == list:
        if len(bestalpha) > 1:
            raise Exception('Returned more than one alpha value, check codes')
        bestalpha = bestalpha[0]
    #import pdb; pdb.set_trace()
    #assert False
    return bestalpha, fit1, size1, alphas


def findalphaD(Ghat, dhat, I, zeroTime, samplerate, numsta, datlenorig, tolerance=None,
               L1=0., L2=0., Tikhratio=[1., 0., 0.]):
    """
    Use discrepancy principle and noise window to find best alpha (tends to find value that
    is too large such that amplitudes are damped)

    Uses Mozorokov discrepancy principle (and bisection method in log-log space?) to find an
    appropriate value of alpha that results in a solution with a fit that is slightly larger than
    the estimated noise level

    Args:
        tolerance (float): how close you want to get to the noise level with the solution
    """
    # Estimate the noise level (use signal before zeroTime)
    dtemp = dhat.copy()[:numsta*datlenorig]  # Trim off any extra zeros
    #lenall = len(dtemp)
    dtemp = np.reshape(dtemp, (numsta, datlenorig))
    if zeroTime is None:
        print('zeroTime not defined, noise estimated from first 100 samples')
        samps = 100
    else:
        samps = int(zeroTime*samplerate)
    temp = dtemp[:, :samps]
    noise = np.sum(np.sqrt(datlenorig * np.std(temp, axis=1)**2))

    # Find ak and bk that yield f(alpha) = ||Gm-d|| - ||noise|| that have f(alpha) values with
    # opposite signs
    templ1 = np.floor(np.log10(np.linalg.norm(Ghat)))
    templ2 = np.arange(templ1-6, templ1)
    ak = 10**templ2[0]
    bk = 10**templ2[-1]
    opposite = False
    fit1 = []
    size1 = []
    alphas = []
    
    Ghat = np.matrix(Ghat)

    Apart = np.dot(Ghat.H, Ghat)
    if type(L2) == float or type(L2) == int:
        L2part = 0.
    else:
        L2part = np.dot(L2.T, L2)
    if type(L1) == float or type(L2) == int:
        L1part = 0.
    else:
        L1part = np.dot(L1.T, L1)

    x = np.squeeze(np.asarray(np.dot(Ghat.H, dhat)))
    
    while opposite is False:
        print(('ak = %s' % (ak,)))
        print(('bk = %s' % (bk,)))
        Aa = Apart+ak**2*(Tikhratio[0]*I + Tikhratio[1]*L1part + Tikhratio[2]*L2part)
        modelak, residuals, rank, s = sp.linalg.lstsq(Aa, x)
        Ab = Apart+bk**2*(Tikhratio[0]*I + Tikhratio[1]*L1part + Tikhratio[2]*L2part)
        modelbk, residuals, rank, s = sp.linalg.lstsq(Ab, x)
        fitak = sp.linalg.norm(np.dot(Ghat, modelak.T)-dhat)
        fitbk = sp.linalg.norm(np.dot(Ghat, modelbk.T)-dhat)
        # Save info on these runs for Lcurve later if desired
        fit1.append(fitak)
        alphas.append(ak)
        size1.append(sp.linalg.norm(Tikhratio[0]*modelak) + sp.linalg.norm(Tikhratio[1]*np.dot(L1part, modelak))\
                     + sp.linalg.norm(Tikhratio[2]*np.dot(L2part, modelak)))
        #sp.linalg.norm(Tikhratio[0]*modelak + Tikhratio[1]*L1part*modelak + Tikhratio[2]*L2part*modelak))#sp.linalg.norm(modelak))
        fit1.append(fitbk)
        alphas.append(bk)
        size1.append(sp.linalg.norm(Tikhratio[0]*modelbk) + sp.linalg.norm(Tikhratio[1]*np.dot(L1part, modelbk))\
                     + sp.linalg.norm(Tikhratio[2]*np.dot(L2part, modelbk)))        
        fak = fitak - noise  # should be negative
        fbk = fitbk - noise  # should be positive
        print(('fak = %s' % (fak,)))
        print(('fbk = %s' % (fbk,)))
        if fak*fbk < 0:
            opposite = True
        if fak > 0:
            ak = 10**(np.log10(ak)-1)
        if fbk < 0:
            bk = 10**(np.log10(bk)+1)

    if tolerance is None:
        tolerance = noise/10.

    # Now use bisection method to find the best alpha value within tolerance
    tol = noise + 100.
    while tol > tolerance:
        # Figure out whether to change ak or bk
        # Compute midpoint (in log units)
        ck = 10**(0.5*(np.log10(ak)+np.log10(bk)))
        Ac = Apart+ck**2*(Tikhratio[0]*I + Tikhratio[1]*L1part + Tikhratio[2]*L2part)
        modelck, residuals, rank, s = sp.linalg.lstsq(Ac, x)
        fitck = sp.linalg.norm(np.dot(Ghat, modelck.T)-dhat)
        fit1.append(fitck)
        alphas.append(ck)
        size1.append(sp.linalg.norm(Tikhratio[0]*modelck) + sp.linalg.norm(Tikhratio[1]*np.dot(L1part, modelck))\
                     + sp.linalg.norm(Tikhratio[2]*np.dot(L2part, modelck)))
        fck = fitck - noise
        print(('ck = %s' % (ck,)))
        print(('fitck = %s' % (fitck,)))
        print(('fck = %s' % (fck,)))
        tol = np.abs(fck)
        if fck*fak < 0:
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
    shp = shp[0]*shp[1]
    dt_temp = np.reshape(dt, shp)
    dtnew_temp = np.reshape(dtnew, shp)
    d_dnew2 = (dt_temp-dtnew_temp)**2
    d2 = dt_temp**2
    VR = (1-(np.sum(d_dnew2)/np.sum(d2))) * 100
    return VR


def back2time(d, df_new, numsta, datlenorig):
    """
    convert data back to the time domain and cut off zero padding
    """
    datlength = int(len(d)/numsta)
    dfrsp = np.reshape(d, (numsta, datlength))
    dfnrsp = np.reshape(df_new, (numsta, datlength))
    #for i in np.arange(len(dfrsp)):
    #    dfrsp[i,:]=symmetric(dfrsp[i,:])
    #    dfnrsp[i,:]=symmetric(dfnrsp[i,:])
    dt = np.real(np.fft.ifft(dfrsp, axis=1))
    dt = dt[0:, 0:datlenorig]
    dtnew = np.real(np.fft.ifft(dfnrsp, axis=1))
    dtnew = dtnew[0:, 0:datlenorig]
    return dt, dtnew


def forward_model(G, model):
    """
    run the forward model (without weights in order to compare to unweighted data)
    """
    #model = sparse.csr_matrix(model)
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
        C = np.zeros((2*len(c)-1, 2*len(c)-1))
        for i in range(2*len(c)-1):
            if i > len(c)-1:
                zros = np.zeros(i+1-len(c))
                p = np.concatenate((zros, cflip, np.zeros(2*len(c))))
            else:
                p = np.concatenate(((cflip[-(i+1):]), np.zeros(2*len(c))))
            p = p[:2*len(c)-1]
            C[i, :] = p.copy()
    else:
        #make it the right size
        C = np.zeros(size)
        for i in range(size[0]):
            if i > len(c)-1:
                zros = np.zeros(i+1-len(c))
                p = np.concatenate((zros, cflip, np.zeros(size[1])))
            else:
                p = np.concatenate(((cflip[-(i+1):]), np.zeros(size[1])))
            p = p[:size[1]]  # cut p to the right size
            C[i, :] = p.copy()
    return C


def datafit(dt, dtnew, st, samplerate, zerotime=0.):
    """
    plot comparision
    """
    tvec = np.arange(0, np.shape(dt)[1]*1/samplerate, 1/samplerate)-zerotime
    fig = plt.figure(figsize=(10, 11))
    offset = 1.*np.amax(np.amax(np.absolute(dt), 0))
    oneline = offset*np.ones((1, np.shape(dt)[1]))
    labels = []
    yticks1 = []
    for i, trace in enumerate(st):
        if i == 0:
            addmat = oneline
        else:
            addmat = np.vstack((addmat, oneline*(i+1)))
        temp = trace.stats.station+'.'+trace.stats.channel+'.'+trace.stats.location+'.'+trace.stats.network
        staname = ('%s - %2.1f km') % (temp, trace.stats.rdist)
        labels.append(staname)
        yticks1.append(-(i+1)*offset)

    ax = fig.add_axes([0.25, 0.05, 0.7, 0.9])
    #.T might be flipping the data upside down...
    ax.plot(np.tile(tvec, (len(st), 1)).T, dt.T-addmat.T, 'k', label='Original')
    ax.plot(np.tile(tvec, (len(st), 1)).T, dtnew.T-addmat.T, 'r', label='Model')
    ax.set_xlim((tvec[0], tvec[-1]))
    ax.set_xlabel('Time (sec)')
    ax.set_yticks(yticks1)
    ax.set_yticklabels(labels)
    ax.set_title('Variance Reduction %2.0f%%' % varred(dt, dtnew))
    redline = mlines.Line2D([], [], color='r', label='Model')
    blackline = mlines.Line2D([], [], color='k', label='Data')
    ax.legend(handles=[blackline, redline], loc='upper right')
    plt.show()
    return fig


def nextpow2(val):
    import math
    temp = math.floor(math.log(val, 2))
    return int(math.pow(2, temp+1))
    pass


def plotinv(Zforce, Nforce, Eforce, tvec, zerotime=0., subplots=False, xlim=None, ylim=None,
            sameY=True, Zupper=None, Zlower=None, Eupper=None, Elower=None, Nupper=None, Nlower=None):
    """
    plot inversion result
    USAGE plotinv(Zforce,Nforce,Eforce,tvec,T0,zerotime=0.,subplots=False,Zupper=None,Zlower=None,Eupper=None,Elower=None,Nupper=None,Nlower=None):
    INPUTS
    [ZEN]force
    tvec =
    T0 = T0 (time delay) used in Green's functions (usually negative)
    zerotime = designated time for event start
    subplots = True, make subplots, False, plot all one one plot
    vline = plot vertical line at t=vline
    [ZEN]upper = upper limit of uncertainties (None if none)
    [ZEN]lower = ditto for lower limit
    OUPUTS
    fig - figure handle
    """
    tvec = tvec - zerotime
    if ylim is None and Zupper is None:
        ylim1 = (np.amin([Zforce.min(), Eforce.min(), Nforce.min()]), np.amax([Zforce.max(),
                         Eforce.max(), Nforce.max()]))
        ylim = (ylim1[0]+0.1*ylim1[0], ylim1[1]+0.1*ylim1[1])  # add 10% on each side to make it look nicer
    elif ylim is None and Zupper is not None:
        ylim1 = (np.amin([Zlower.min(), Elower.min(), Nlower.min()]), np.amax([Zupper.max(),
                         Eupper.max(), Nupper.max()]))
        ylim = (ylim1[0]+0.1*ylim1[0], ylim1[1]+0.1*ylim1[1])  # add 10% on each side to make it look nicer
    if subplots:
        fig = plt.figure(figsize=(14, 9))
        ax1 = fig.add_subplot(311)
        ax1.plot(tvec, Zforce, 'b')
        ax1.grid(True)
        ax1.set_title('Up')
        ax2 = fig.add_subplot(312)  # ,sharex=ax1)
        ax2.plot(tvec, Nforce, 'r')
        ax2.grid(True)
        ax2.set_title('North')
        ax3 = fig.add_subplot(313)  # ,sharex=ax1)
        ax3.plot(tvec, Eforce, 'g')
        ax3.grid(True)
        ax3.set_title('East')
        if sameY or ylim is not None:
            ax1.set_ylim(ylim)
            ax2.set_ylim(ylim)
            ax3.set_ylim(ylim)
        x = np.concatenate((tvec, tvec[::-1]))
        if Zupper is not None and Zlower is not None:
            y = np.concatenate((Zlower, Zupper[::-1]))
            poly = plt.Polygon(list(zip(x, y)), facecolor='b', edgecolor='none', alpha=0.2)
            ax1.add_patch(poly)
        if Nupper is not None and Nlower is not None:
            y = np.concatenate((Nlower, Nupper[::-1]))
            poly = plt.Polygon(list(zip(x, y)), facecolor='r', edgecolor='none', alpha=0.2)
            ax2.add_patch(poly)
        if Eupper is not None and Elower is not None:
            y = np.concatenate((Elower, Eupper[::-1]))
            poly = plt.Polygon(list(zip(x, y)), facecolor='g', edgecolor='none', alpha=0.2)
            ax3.add_patch(poly)
        if xlim:
            ax1.set_xlim(xlim)
            ax2.set_xlim(xlim)
            ax3.set_xlim(xlim)
        ax2.set_ylabel('Force (N)')
    else:
        fig = plt.figure(figsize=(14, 4))
        ax = fig.add_subplot(111)
        ax.plot(tvec, Zforce, 'b', label='Up')
        ax.plot(tvec, Nforce, 'r', label='North')
        ax.plot(tvec, Eforce, 'g', label='East')
        x = np.concatenate((tvec, tvec[::-1]))
        if Zupper is not None and Zlower is not None:
            y = np.concatenate((Zlower, Zupper[::-1]))
            poly = plt.Polygon(list(zip(x, y)), facecolor='b', edgecolor='none', alpha=0.2)
            ax.add_patch(poly)
        if Nupper is not None and Nlower is not None:
            y = np.concatenate((Nlower, Nupper[::-1]))
            poly = plt.Polygon(list(zip(x, y)), facecolor='r', edgecolor='none', alpha=0.2)
            ax.add_patch(poly)
        if Eupper is not None and Elower is not None:
            y = np.concatenate((Elower, Eupper[::-1]))
            poly = plt.Polygon(list(zip(x, y)), facecolor='g', edgecolor='none', alpha=0.2)
            ax.add_patch(poly)
        if xlim:
            ax.set_xlim(xlim)
        ax.legend(loc='upper right')
        ax.grid(True)
        ax.set_ylabel('Force (N)')
        ax.set_ylim(ylim)
    plt.xlabel('Time (sec')
    plt.show()
    return fig


def plotangmag(Zforce, Nforce, Eforce, tvec, zerotime=0., subplots=False, xlim=None, ylim=None,
               sameY=True, Zupper=None, Zlower=None, Eupper=None, Elower=None, Nupper=None,
               Nlower=None):
    """
    plot angles and magnitudes of inversion result
    USAGE plotinv(Zforce,Nforce,Eforce,tvec,T0,zerotime=0.,subplots=False,Zupper=None,Zlower=None,Eupper=None,Elower=None,Nupper=None,Nlower=None):
    INPUTS
    [ZEN]force
    tvec =
    T0 = T0 (time delay) used in Green's functions (usually negative)
    zerotime = designated time for event start
    subplots = True, make subplots, False, plot all one one plot
    vline = plot vertical line at t=vline
    [ZEN]upper = upper limit of uncertainties (None if none)
    [ZEN]lower = ditto for lower limit
    OUPUTS
    fig - figure handle
    """
    tvec = tvec - zerotime
    if ylim is None and Zupper is None:
        ylim1 = (np.amin([Zforce.min(), Eforce.min(), Nforce.min()]), np.amax([Zforce.max(), Eforce.max(), Nforce.max()]))
        ylim = (ylim1[0]+0.1*ylim1[0], ylim1[1]+0.1*ylim1[1])  # add 10% on each side to make it look nicer
    elif ylim is None and Zupper is not None:
        ylim1 = (np.amin([Zlower.min(), Elower.min(), Nlower.min()]), np.amax([Zupper.max(), Eupper.max(), Nupper.max()]))
        ylim = (ylim1[0]+0.1*ylim1[0], ylim1[1]+0.1*ylim1[1])  # add 10% on each side to make it look nicer
        fig = plt.figure(figsize=(10, 10))

    # Plot the inversion result in the first one
    ax = fig.add_subplot(411)
    ax.plot(tvec, Zforce, 'b', label='Up')
    ax.plot(tvec, Nforce, 'r', label='North')
    ax.plot(tvec, Eforce, 'g', label='East')
    x = np.concatenate((tvec, tvec[::-1]))
    if Zupper is not None and Zlower is not None:
        y = np.concatenate((Zlower, Zupper[::-1]))
        poly = plt.Polygon(list(zip(x, y)), facecolor='b', edgecolor='none', alpha=0.2)
        ax.add_patch(poly)
    if Nupper is not None and Nlower is not None:
        y = np.concatenate((Nlower, Nupper[::-1]))
        poly = plt.Polygon(list(zip(x, y)), facecolor='r', edgecolor='none', alpha=0.2)
        ax.add_patch(poly)
    if Eupper is not None and Elower is not None:
        y = np.concatenate((Elower, Eupper[::-1]))
        poly = plt.Polygon(list(zip(x, y)), facecolor='g', edgecolor='none', alpha=0.2)
        ax.add_patch(poly)
    if xlim:
        ax.set_xlim(xlim)
    ax.legend(loc='upper right')
    ax.grid(True)
    ax.set_ylabel('Force (N)')
    ax.set_ylim(ylim)

    # Plot the magnitudes in second one
    ax1 = fig.add_subplot(412)
    Mag = np.linalg.norm(list(zip(Zforce, Eforce, Nforce)), axis=1)
    MagU = np.linalg.norm(list(zip(np.maximum(np.abs(Zupper), np.abs(Zlower)),
                                   np.maximum(np.abs(Eupper), np.abs(Elower)),
                                   np.maximum(np.abs(Nupper), np.abs(Nlower)))), axis=1)
    MagL = np.linalg.norm(list(zip(np.minimum(np.abs(Zupper), np.abs(Zlower)),
                                   np.minimum(np.abs(Eupper), np.abs(Elower)),
                                   np.minimum(np.abs(Nupper), np.abs(Nlower)))), axis=1)
    ax1.plot(tvec, Mag, 'k', label='best')
    ax1.plot(tvec, MagL, 'r', label='lower')
    ax1.plot(tvec, MagU, 'r', label='upper')
    #ax1.legend(loc='upper right')
    ax1.grid(True)
    ax1.set_ylabel('Force (N)')
    if xlim:
        ax1.set_xlim(xlim)

    # Plot the horizontal azimuth
    ax2 = fig.add_subplot(413)
    tempang = (180/np.pi)*np.arctan2(Nforce, Eforce)-90  # get angle counterclockwise relative to N
    #any negative values, add 360
    for i, temp in enumerate(tempang):
        if temp < 0:
            tempang[i] = temp+360
    tempangU = (180/np.pi)*np.arctan2(Nupper, Eupper)-90
    for i, temp in enumerate(tempangU):
        if temp < 0:
            tempangU[i] = temp+360
    tempangL = (180/np.pi)*np.arctan2(Nlower, Elower)-90
    for i, temp in enumerate(tempangL):
        if temp < 0:
            tempangL[i] = temp+360
    # now flip to clockwise to get azimuth
    Haz = 360-tempang
    #HazU = 360-tempangU
    #HazL = 360-tempangL
    ax2.plot(tvec, Haz)
    #ax2.plot(tvec, HazU, 'r')
    #ax2.plot(tvec, HazL, 'r')
    ax2.grid(True)
    ax2.set_ylabel('Azimuth (deg CW from N)')
    if xlim:
        ax2.set_xlim(xlim)

    #Plot the vertical angle
    ax3 = fig.add_subplot(414)
    Vang = (180/np.pi)*np.arctan(Zforce/np.sqrt(Nforce**2+Eforce**2))
    #VangU = (180/np.pi)*np.arctan(Zlower/np.sqrt(Nupper**2+Elower**2))
    #VangL = (180/np.pi)*np.arctan(Zupper/np.sqrt(Nlower**2+Eupper**2))
    ax3.plot(tvec, Vang)
    ax3.grid(True)
    ax3.set_ylabel('Vertical angle (deg)')
    #ax3.plot(tvec, VangU, 'r')
    #ax3.plot(tvec, VangL, 'r')
    if xlim:
        ax3.set_xlim(xlim)
    plt.xlabel('Time (sec')
    plt.show()

    return fig, Mag, MagU, MagL, Vang, Haz


def curvature(x, y, negslope=True):
    """
    Estimate the radius of curvature for each point on line to find corner of L-curve
    
    Args:
        x (array): x points
        y (array): y points

    Returns:
        radius of curvature for each point (ends will be nan)
    """

    #FOR EACH SET OF THREE POINTS, FIND RADIUS OF CIRCLE THAT FITS THEM - IGNORE ENDS
    R_2 = np.ones(len(x))*float('inf')  # end numbers should be infinity because is straight line
    for i in range(1, len(R_2)-1):
        xsub = x[i-1:i+2]
        ysub = y[i-1:i+2]
        m1 = -1/((ysub[0]-ysub[1])/(xsub[0]-xsub[1]))  # slope of bisector of first segment
        m2 = -1/((ysub[1]-ysub[2])/(xsub[1]-xsub[2]))  # slope of bisector of second segment
        b1 = ((ysub[0]+ysub[1])/2)-m1*((xsub[0]+xsub[1])/2)  # compute b for first bisector
        b2 = ((ysub[1]+ysub[2])/2)-m2*((xsub[1]+xsub[2])/2)  # compute b for second bisector

        Xc = (b1-b2)/(m2-m1)  # find intercept point of bisectors
        Yc = b2 + m2*Xc

        R_2[i] = np.sqrt((xsub[0]-Xc)**2 + (ysub[0]-Yc)**2)  # get distance from any point to intercept of bisectors to get radius

    return R_2


def rotate2ZNE(data_1, azimuth_1, dip_1, data_2, azimuth_2, dip_2, data_3, azimuth_3, dip_3):
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
    dip_rotation_matrix = np.cos(dip) * \
        np.matrix(((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))) + \
        (1 - np.cos(dip)) * np.matrix(((c1 * c1, c1 * c2, c1 * c3),
                                      (c2 * c1, c2 * c2, c2 * c3),
                                      (c3 * c1, c3 * c2, c3 * c3))) + \
        np.sin(dip) * np.matrix(((0, -c3, c2), (c3, 0, -c1), (-c2, c1, 0)))
    # Do the same for the azimuth.
    c1 = -1.0
    c2 = 0.0
    c3 = 0.0
    azimuth_rotation_matrix = np.cos(azimuth) * \
        np.matrix(((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))) + \
        (1 - np.cos(azimuth)) * np.matrix(((c1 * c1, c1 * c2, c1 * c3),
                                          (c2 * c1, c2 * c2, c2 * c3),
                                          (c3 * c1, c3 * c2, c3 * c3))) + \
        np.sin(azimuth) * np.matrix(((0, -c3, c2), (c3, 0, -c1), (-c2, c1, 0)))
    # Now simply rotate a north pointing unit vector with both matrixes.
    temp = np.dot(azimuth_rotation_matrix, [[0.0], [-1.0], [0.0]])
    return np.array(np.dot(dip_rotation_matrix, temp)).ravel()


def saverun(filename, model, Zforce, Eforce, Nforce, tvec, dt, dtnew, stproc, zeroTime,
            ZforceU=None, ZforceL=None, EforceU=None, EforceL=None, NforceU=None, NforceL=None):
    """
    Save the results of an inversion, regular or jackknife
    Makes a dictionary and saves that as a pickle
    USAGE: saverun(filename, model, Zforce, Eforce, Nforce, tvec, dt, dtnew, stproc, zeroTime,
                   ZforceU=None, ZforceL=None, EforceU=None, EforceL=None, NforceU=None, NforceL=None)
    RETURN nothing but saves file as filename
    """
    f = open(filename, 'wb')
    result = {'model': model, 'tvec': tvec, 'Zforce': Zforce, 'Eforce': Eforce, 'Nforce': Nforce,
              'ZforceU': ZforceU, 'ZforceL': ZforceL, 'EforceU': EforceU, 'EforceL': EforceL,
              'NforceU': NforceU, 'NforceL': NforceL, 'dt': dt, 'dtnew': dtnew, 'stproc': stproc,
              'zeroTime': zeroTime, }
    pickle.dump(result, f)
    f.close()
    return


def readrun(filename):
    """
    Read the results of an inversion, regular or jackknife, and outputs them as original variable 
        names instead of in dictionary
    USAGE: readrun(filename)
    RETURN: model, tvec, Zforce, Eforce, Nforce, ZforceU, ZforceL, EforceU, EforceL, NforceU,
        NforceL, dt, dtnew, stproc, zeroTime
    """
    f = open(filename, 'rb')
    result = pickle.load(f)
    f.close()
    model = result['model']
    tvec = result['tvec']
    Zforce = result['Zforce']
    Eforce = result['Eforce']
    Nforce = result['Nforce']
    ZforceU = result['ZforceU']
    ZforceL = result['ZforceL']
    EforceU = result['EforceU']
    EforceL = result['EforceL']
    NforceU = result['NforceU']
    NforceL = result['NforceL']
    dt = result['dt']
    dtnew = result['dtnew']
    stproc = result['stproc']
    zeroTime = result['zeroTime']
    return(model, tvec, Zforce, Eforce, Nforce, ZforceU, ZforceL, EforceU, EforceL, NforceU,
           NforceL, dt, dtnew, stproc, zeroTime)

#def symmetric(data):
#    """
#    make a vector symmetric before doing ifft so it turns out real
#    """
#    datlen = len(data)
#    if datlen%2 == 1: #odd
#        out = np.hstack((0.,data[1:datlen/2],np.conjugate(data[1:datlen/2])[::-1]))
#    else: #even
#        out = np.hstack((0.,data[1:datlen/2],0.,np.conjugate(data[1:datlen/2])[::-1]))
#    if datlen != len(out):
#        print 'input data different length than output - fix something!'
#        return
#    else:
#        return data
