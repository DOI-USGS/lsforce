#!/usr/bin/env python

"""
Functions for outputting info from the database
"""

import findsta
from reviewData import reviewData
import numpy as np
from obspy import Stream, UTCDateTime
from obspy.clients.fdsn import Client as FDSN_Client
import glob
import urllib2
import os

# bufferperc=0.1
# HFlims=(1., 5.)
# HFoutput='VEL'
# LPlims=(20., 60.)
# LPoutput='DISP'
# mintraces=5
# maxtraces=15
# path='/Users/kallstadt/LSseis/landslideDatabase'
# database='/Users/kallstadt/LSseis/landslideDatabase/lsseis.db'


def summary_file():
    """
    Output a csv file summarizing the database (event table primarily)
    """
    pass


def grab_data(event_id, bufferperc=0.1, numstas=5, path='/Users/kallstadt/LSseis/landslideDatabase',
              detectHF=True, detectLP=True, minradius=0., maxradius=None, Zonly=False,
              database='/Users/kallstadt/LSseis/landslideDatabase/lsseis.db'):

    """
    Pulls data from IRIS and other sources and links station correction info
    :param event_id: Integer specifying which event to review
    :param bufferperc: Buffer percentage of duration (end time - start time) (makes viewing signal easier)
    :param numstas: Number of stations to show (shows numstas closest stations, all channels)
    :param database: Full file path of database file location
    :param path: Path to location of sac files (upstream from relative file paths listed in database)

    """
    # Get event info
    evDict = findsta.getEventInfo(event_id, database=database)
    print(('Now grabbing requested data for Eid %s - %s') % (event_id, evDict['Name']))

    buffer_sec = (evDict['EndTime'] - evDict['StartTime']) * bufferperc

    # Find stations where detect_HF=1 or detect_LP=1 for this event
    statemp = findsta.getStaInfo(event_id, database=database, minradius=minradius, maxradius=maxradius)
    # Loop through and keep only the ones that were detected on HF and/or LP
    keep = []
    for k, v in statemp.items():
        if detectHF:
            if v['detect_HF'] == 1:
                keep.append(v)
            continue
        else:
            keep.append(v)
            continue
        if detectLP:
            if v['detect_LP'] == 1:
                keep.append(v)
        else:
            keep.append(v)

    # Then keep only the numchans closest stations (all channels)
    indx, dists = zip(*[[i, k['stasource_radius_km']] for i, k in enumerate(keep)])
    dists = reviewData.unique_list(sorted(dists))
    # Take first x unique distances (numchans)
    unique_dist = reviewData.unique_list(dists)
    these = unique_dist[:numstas]

    # Take those ones from the keep dictionary
    staDict = []
    for k in keep:
        if k['stasource_radius_km'] in these:
            staDict.append(k)

    datlocs = reviewData.unique_list([k['source'] for k in staDict])

    sttemp = Stream()
    if 'sac' in evDict['DatLocation']:
        if Zonly:
            chanuse = 'BHZ,EHZ,HHZ'
        else:
            chanuse = '*'
        datloc1 = evDict['DatLocation'].split(',')
        datloc1 = [x.strip() for x in datloc1 if 'sac' in x]
        for datl in datloc1:
            fullpath = os.path.join(path, datl.split(':')[1])
            filenames = glob.glob(fullpath)
            if len(filenames) > 0:
                if 'AlaskaData' in datl:  # Need to do some cheating to attach response info if Iliamna sac data - attach oldest response info available at IRIS for each station
                    # Only keep filenames of stations we want to read in
                    newfilenames = []
                    namelist = [k['Name'] for k in staDict]
                    for filen in filenames:
                        nam = filen.split('/')[-1].split('.')[0]
                        if nam in namelist:
                            newfilenames.append(filen)
                    stsac = reviewData.getdata_sac(newfilenames, attach_response=True, starttime=evDict['StartTime']-buffer_sec,
                                                   endtime=evDict['EndTime']+buffer_sec, chanuse=chanuse)
                    stsac = Stream([trace for trace in stsac if trace.max() != 0.0])  # Get rid of any empty ones
                    originalstt = []
                    # Get earliest start time from IRIS for AK and AV
                    url1 = ('http://service.iris.edu/fdsnws/station/1/query?network=AV&level=station&format=text&nodata=404')
                    url2 = ('http://service.iris.edu/fdsnws/station/1/query?network=AK&level=station&format=text&nodata=404')
                    f = urllib2.urlopen(url1)
                    file1 = f.read()
                    lines1 = [line.split('|') for line in file1.split('\n')[1:]]
                    stdict = {}
                    f.close()
                    f = urllib2.urlopen(url2)
                    file2 = f.read()
                    lines2 = [line.split('|') for line in file2.split('\n')[1:]]
                    f.close()
                    for line in lines1[:-1]:
                        stdict[line[1]] = (line[0], UTCDateTime(line[6]))
                    for line in lines2[:-1]:
                        if line[1] in stdict.keys():
                            # Take the minimum if shows up for both network codes
                            stdict[line[1]] = (line[0], np.min([stdict[line[1]], UTCDateTime(line[6])]))
                        else:
                            stdict[line[1]] = (line[0], UTCDateTime(line[6]))

                    for temp in stsac:
                        originalstt.append(temp.stats.starttime)
                        if temp.stats.channel[0] == 'S':
                            temp.stats.channel = 'E'+temp.stats.channel[1:]
                        # Make sure network code is consistent with IRIS
                        if temp.stats.station in stdict:
                            temp.stats.network = stdict[temp.stats.station][0]
                        else:
                            print 'could not attach response info for %s, station correction will not work' % temp.stats.station
                            continue
                        if 'response' not in temp.stats.keys():
                            # Change starttime of temp to earliest start time at IRIS for this station
                            try:
                                temp.stats.starttime = stdict[temp.stats.station][1] + 86400.
                            except:
                                print 'could not attach response info for %s, station correction will not work' % temp.stats.station
                    # Get responses from IRIS
                    client = FDSN_Client('IRIS')
                    akstas = ','.join([tr.stats.station for tr in stsac])
                    inv1 = client.get_stations(network='AK,AV', station=akstas, level="response")
                    stsac.attach_response(inv1)
                    # Change back start times
                    for i, trace in enumerate(stsac):
                        trace.stats.starttime = originalstt[i]
                else:
                    stsac = reviewData.getdata_sac(filenames, attach_response=True, starttime=evDict['StartTime']-buffer_sec,
                                                   endtime=evDict['EndTime']+buffer_sec, chanuse=chanuse)
                sttemp += stsac

                for trace in sttemp:
                    if '--' in trace.stats.location:
                        trace.stats.location = ''

    # remove any from staDict that were already loaded before moving on to save time and avoid duplicates
    staDict2 = [sta for sta in staDict if sta['Name'] not in [trace.stats.station for trace in sttemp]]

    if 'IRIS' in evDict['DatLocation'] or 'IRIS' in datlocs:
        stalist = []
        for k in staDict2:
            if Zonly is True:
                if k['Channel'][-1].lower() != 'z':
                    continue
            if 'IRIS' in k['source']:
                if 'Iliamna' not in evDict['DatLocation']:
                    stalist.append((k['Name'], k['Channel'], k['Network'], '*'))
                else:
                    if 'AV' not in k['Network'] and 'AK' not in k['Network']:
                        stalist.append((k['Name'], k['Channel'], k['Network'], '*'))
        if len(stalist) != 0:
            try:
                sttemp += reviewData.getdata_exact(stalist, evDict['StartTime'] - buffer_sec, evDict['EndTime'] + buffer_sec,
                                                   attach_response=True, clientname='IRIS')
            except Exception as e:
                print(e)
    if 'NCEDC' in evDict['DatLocation'] or 'NCEDC' in datlocs:
        stalist = [(k['Name'], k['Channel'], k['Network'], '*') for k in staDict2 if 'NCEDC' in k['source']]
        if len(stalist) != 0:
            sttemp += reviewData.getdata_exact(stalist, evDict['StartTime'] - buffer_sec, evDict['EndTime'] + buffer_sec,
                                               attach_response=True, clientname='NCEDC')

    # Delete any that are not in list
    nets, stas, chans = zip(*[[k['Network'], k['Name'], k['Channel']] for k in staDict])  # Exclude location code because inconsistencies here can cause data to not be found
    st = Stream()
    for n1, s1, c1 in zip(nets, stas, chans):
        st += sttemp.select(station=s1, network=n1, channel=c1)
    # Attach distaz
    st = findsta.attach_distaz(st, evDict['Latitude'], evDict['Longitude'], database=database)
    st = st.sort(keys=['rdist', 'channel'])

    # Preprocess
    st.detrend('demean')

    return st


def make_figures(st, bufferperc=None, raw=True, HF=True, LP=True, timeseries=True, spectrograms=True, spectra=True,
                 HFlims=(1., 5.), HFoutput='VEL', LPlims=(20., 60.), LPoutput='DISP'):
    """
    Make static figures for each event
    :param HFlims: tuple or list of lower and upper frequency limits in Hz for HF filtering (1-5 Hz standard)
    :param HFoutput: Output type for HF station correction (VEL standard)
    UPDATE THESE:
    spectra = multitaper, equiv to PSD
    """
    # build cosine filter that will be used
    cosfiltHF = (0.5*HFlims[0], HFlims[0], HFlims[1], 2*HFlims[1])
    cosfiltLP = (0.5/LPlims[1], 1/LPlims[1], 1/LPlims[0], 2/LPlims[0])

    if bufferperc is None:
        taper = 0.05
    else:
        taper = bufferperc/4.

    # Do corrections for LP, HF
    if LP:
        stLP = st.copy()
        stLP.taper(max_percentage=taper, type='cosine')
        try:
            stLP.remove_response(output=LPoutput, pre_filt=cosfiltLP)
        except:
            temp = Stream()
            for i, trace in enumerate(stLP):
                try:
                    trace.remove_response(output=LPoutput, pre_filt=cosfiltLP)
                    temp = temp + trace
                except:
                    print 'Failed to remove response for %s, deleting this station' % (trace.stats.station + trace.stats.channel,)
            stLP = temp.copy()

    if HF:
        stHF = st.copy()
        stHF.taper(max_percentage=taper, type='cosine')
        try:
            stHF.remove_response(output=HFoutput, pre_filt=cosfiltHF)
        except:
            temp = Stream()
            for i, trace in enumerate(stHF):
                try:
                    trace.remove_response(output=HFoutput, pre_filt=cosfiltHF)
                    temp = temp + trace
                except:
                    print 'Failed to remove response for %s, deleting this station' % (trace.stats.station + trace.stats.channel,)
            stHF = temp.copy()

    # Figure out the colors, make sure same station has same color
    cycle = ['r', 'b', 'g', 'm', 'k']
    colors1 = []
    stas = reviewData.unique_list([trace.stats.station for trace in st])
    for trace in st:
        ind = stas.index(trace.stats.station)
        colors1 += cycle[ind]

    # Now make figures
    if timeseries:
        figraw = reviewData.recsec(st, norm=True, maxtraces=15, quickdraw=True, figsize=(13, 14), colors=colors1,
                                   labelsize=14, addscale=False, unitlabel=None, convert=1., labelquickdraw=False)
        if HF:
            figHF = reviewData.recsec(stHF, norm=False, maxtraces=15, quickdraw=True, figsize=(13, 14), colors=colors1,
                                      labelsize=14, addscale=True, unitlabel='m/s', convert=1., labelquickdraw=False)
        else:
            figHF = None
        if LP:
            figLP = reviewData.recsec(stLP, norm=False, maxtraces=15, quickdraw=True, figsize=(13, 14), colors=colors1,
                                      labelsize=14, addscale=True, unitlabel='m', convert=1., labelquickdraw=False)
        else:
            figLP = None
    else:
        figraw = None
        figHF = None
        figLP = None

    # Make spectrograms and spectra
    if spectrograms or spectra:
        stspec = st.copy()
        try:
            stspec.remove_sensitivity()
        except:
            temp = Stream()
            for i, trace in enumerate(stspec):
                try:
                    trace.remove_sensitivity()
                    temp = temp + trace
                except:
                    print 'Failed to remove sensitivity for %s, deleting this station' % (trace.stats.station + trace.stats.channel,)

    if spectrograms:
        specg = reviewData.make_spectrogram(stspec, log1=True, maxtraces=15, maxPower=None, minPower=None, freqmax=None, labelsize=14)
    else:
        specg = None
    if spectra:
        mt = reviewData.make_multitaper(st, number_of_tapers=None, time_bandwidth=4., sine=False, recsec=True, colors1=colors1,
                                        logx=True, logy=True, xunits='Hz', xlim=None, yunits='$m^2/s$', ylim=(10.**-3, 10.**9))
    else:
        mt = None

    return figraw, figHF, figLP, specg, mt
