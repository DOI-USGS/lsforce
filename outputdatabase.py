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
from cleandatabase import get_eids
import sqlite3 as lite
import datetime
import matplotlib.pyplot as plt

# bufferperc=0.1
# HFlims=(1., 5.)
# HFoutput='VEL'
# LPlims=(20., 60.)
# LPoutput='DISP'
# mintraces=5
# maxtraces=15
# path='/Users/kallstadt/LSseis/landslideDatabase'
# database='/Users/kallstadt/LSseis/landslideDatabase/lsseis.db'


def waveformfig_db(eids=None, numstas=5, bufferperc=0.15, database='/Users/kallstadt/LSseis/landslideDatabase/lsseis.db',
                   placefigs='/Users/kallstadt/LSseis/landslideDatabase/InfoFiles/waveformfigs', relplacefigs='InfoFiles/waveformfigs',
                   savedat=False, folderdat='data', reloadfile=False, Zonly=False, addtodb=False, loadfromfile=True,
                   minradius=0., maxradius=None, removeoutliers=False, raw=True, HF=True, LP=True, timeseries=True,
                   spectrograms=True, spectra=True, merge=True, pad=True, fill_value=0., detrend='demean'):
    """
    Make standard figures for each event and put them in the database
    """

    try:
        os.mkdir(placefigs)
    except Exception as e:
        print(e)

    if eids is None:
        eids = get_eids(database)

    if type(eids) is int:
        eids = [eids]

    for e1 in eids:
        evDict = findsta.getEventInfo(e1)
        detectLP = False
        if evDict['LPpotential'] == 1:
            detectLP = True
        st = grab_data(e1, both=True, database=database, savedat=savedat, folderdat=folderdat, minradius=minradius, maxradius=maxradius,
                       bufferperc=bufferperc, numstas=numstas, Zonly=Zonly, reloadfile=reloadfile, loadfromfile=loadfromfile,
                       merge=merge, pad=pad, fill_value=fill_value, detrend=detrend)
        if len(st) == 0:
            print('No data for event %d, continuing to next eid' % e1)
            continue
        print('got data')
        # remove duplicates, take higher sample rate station

        unilist = reviewData.unique_list([trace.stats.station for trace in st])

        for sta in unilist:
            shortlist = st.select(station=sta)
            if len(shortlist) > 3:
                samprates = [tr.stats.sampling_rate for tr in shortlist]
                minrate = np.min(samprates)
                if np.min(samprates) != np.max(samprates):
                    for tr in st.select(station=sta, sampling_rate=minrate):
                        if tr.stats.station == sta and tr.stats.sampling_rate == minrate:
                            st.remove(tr)
                # loccodes = np.sort(reviewData.unique_list([trace.stats.location for trace in shortlist]))
                # if len(loccodes) > 1:
                #     for loc in loccodes:
                #         if loc == '' or loc == '00':
                #             pass
                #         else:
                #             for tr in st.select(station=sta, location=loc):
                #                 st.remove(tr)
        if HF:
            # Make high frequency figs (without taper)
            figraw, figHF, junk, specg, junk2 = make_figures(st, bufferperc=None, raw=True, HF=True, LP=False, spectrograms=True, spectra=False,
                                                             taper=None, detrend=detrend, startline=evDict['StartTime'], endline=evDict['EndTime'],
                                                             removeoutliers=removeoutliers)
        if raw and not HF:
            figraw, figHF, junk, specg, junk2 = make_figures(st, bufferperc=None, raw=True, HF=False, LP=False, spectrograms=False, spectra=False,
                                                             taper=None, detrend=detrend, startline=evDict['StartTime'], endline=evDict['EndTime'],
                                                             removeoutliers=removeoutliers)
            figHF = None

        if not spectrograms:
            specg = None

        print('made high frequency plots')
        # Make lp figs and spectra (with taper)
        if detectLP and LP:
            stLP = grab_data(e1, detectHF=None, detectLP=True, database=database, savedat=savedat, folderdat=folderdat + 'LP', minradius=minradius, maxradius=maxradius,
                             bufferperc=bufferperc, numstas=numstas, Zonly=Zonly, reloadfile=reloadfile, loadfromfile=loadfromfile,
                             merge=merge, pad=pad, fill_value=fill_value, detrend=detrend)
            for tr in stLP.select(channel='BDF'):
                stLP.remove(tr)
            unilist = reviewData.unique_list([trace.stats.station for trace in stLP])
            for sta in unilist:
                shortlist = stLP.select(station=sta)
                if len(shortlist) > 3:
                    samprates = [tr.stats.sampling_rate for tr in shortlist]
                    minrate = np.min(samprates)
                    if np.min(samprates) != np.max(samprates):
                        for tr in stLP.select(station=sta, sampling_rate=minrate):
                            if tr.stats.station == sta and tr.stats.sampling_rate == minrate:
                                stLP.remove(tr)

            if LP and timeseries and spectra:
                # make displacement spectra too...
                junk1, junk2, figLP, junk3, mtLP = make_figures(stLP, bufferperc=bufferperc, raw=False, HF=False, LP=True, spectrograms=False, spectra=True,
                                                                taper='bufferperc', detrend='linear', startline=evDict['StartTime'], endline=evDict['EndTime'],
                                                                speccorr='DISP', removeoutliers=removeoutliers)
            elif LP and not timeseries and spectra:
                junk1, junk2, figLP, junk3, mtLP = make_figures(stLP, bufferperc=bufferperc, raw=False, HF=False, LP=False, spectrograms=False, spectra=True,
                                                                taper='bufferperc', detrend='linear', startline=evDict['StartTime'], endline=evDict['EndTime'],
                                                                speccorr='DISP', removeoutliers=removeoutliers)
            else:
                figLP = None
                mtLP = None
        else:
            figLP = None
            mtLP = None

        if spectra and (HF or LP):
            if detrend is None:
                junk1, junk2, junk, junk3, mt = make_figures(st, bufferperc=bufferperc, raw=False, HF=False, LP=False, spectrograms=False, spectra=True,
                                                             taper='bufferperc', detrend=None, removeoutliers=removeoutliers)
            else:
                junk1, junk2, junk, junk3, mt = make_figures(st, bufferperc=bufferperc, raw=False, HF=False, LP=False, spectrograms=False, spectra=True,
                                                             taper='bufferperc', detrend='linear', removeoutliers=removeoutliers)
        else:
            mt = None

        print('finished the rest of the plots')
        # make file name prefix
        evname = evDict['Name'].strip().replace(' ', '').replace(',', '')
        evdat = evDict['StartTime'].strftime('%d%b%Y')
        temp = '%s_%s%s' % (e1, evname, evdat)

        # Save them all

        if figraw is not None:
            newfileR = '%s_waveforms_raw.png' % temp
            figraw.savefig(os.path.join(placefigs, newfileR))
            print('saved raw file')

        if figHF is not None:
            newfileHF = '%s_waveforms_1-5Hz.png' % temp
            figHF.savefig(os.path.join(placefigs, newfileHF))
            print('saved HF file')

        if figLP is not None:
            newfileLP = '%s_waveforms_20-60sec.png' % temp
            figLP.savefig(os.path.join(placefigs, newfileLP))
            print('saved LP file')

        if mtLP is not None:
            newfileSGLP = '%s_displspectra.png' % temp
            mtLP.savefig(os.path.join(placefigs, newfileSGLP))
            print('saved LP spec file')

        if specg is not None:
            newfileSG = '%s_spectrogram.png' % temp
            specg.savefig(os.path.join(placefigs, newfileSG))
            print('saved specg')

        if mt is not None:
            newfileSP = '%s_spectra.png' % temp
            mt.savefig(os.path.join(placefigs, newfileSP))
            print('saved spectra')

        if addtodb:
            # Connect to database
            connection = None
            connection = lite.connect(database)
            cursor = connection.cursor()

            # Get max photo id currently existing
            idlist = cursor.execute("""SELECT photo_id FROM photos""")
            phoids = idlist.fetchall()
            phoids = [id1[0] for id1 in phoids]
            currentpid = np.max(phoids) + 1
            date = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')

            # Put in database
            with connection:
                if figraw is not None:
                    connection.execute('INSERT INTO photos(photo_id, event_id, photographer, description, date, file_extension) VALUES (?,?,?,?,?,?)',
                                       (currentpid, e1, 'Kate Allstadt', 'Raw waveforms from closest stations. Data was demeaned but no other processing was applied. Dashed lines show approximate start and end times of event (from database). Figure generation is automated and there may be artifacts present.', date, os.path.join(relplacefigs, newfileR)))
                    currentpid += 1
                if figHF is not None:
                    connection.execute('INSERT INTO photos(photo_id, event_id, photographer, description, date, file_extension) VALUES (?,?,?,?,?,?)',
                                       (currentpid, e1, 'Kate Allstadt', 'High frequency (1-5 Hz) ground velocity waveforms for station response from closest stations with detections. Dashed lines show approximate start and end times of event (from database). Waveforms are corrected for station response after removing the mean. Figure generation is automated and there may be artifacts present. Clipped waveforms have not been removed before station correction so check raw plots before interpreting high amplitude signals.', date, os.path.join(relplacefigs, newfileHF)))
                    currentpid += 1
                if figLP is not None:
                    connection.execute('INSERT INTO photos(photo_id, event_id, photographer, description, date, file_extension) VALUES (?,?,?,?,?,?)',
                                       (currentpid, e1, 'Kate Allstadt', 'Long period (20-60 sec) ground displacement waveforms from closest stations with detections. Dashed lines show approximate start and end times of event (from database, based on high frequency waveform). Waveforms are corrected for station response after detrending using a linear fit and tapering. Zero mean cosine filter applied with station correction may cause acausal arrivals. Figure generation is automated and there may be artifacts present. Signals with clipped waveforms and data gaps have not been removed so check raw figures before interpreting.', date, os.path.join(relplacefigs, newfileLP)))
                    currentpid += 1
                if mtLP is not None:
                    connection.execute('INSERT INTO photos(photo_id, event_id, photographer, description, date, file_extension) VALUES (?,?,?,?,?,?)',
                                       (currentpid, e1, 'Kate Allstadt', 'Displacement spectrum from closest broadband stations with detections, corrected for station response with no prefiltering after detrending using a linear fit and tapering. Note that long period signals often contain noise unrelated to the event signal of interest, especially on horizontal components (e.g. tilt), so interpret with care. Figure generation is automated and there may be artifacts present. Signals with clipped waveforms and data gaps have not been removed so check raw figures before interpreting.', date, os.path.join(relplacefigs, newfileSGLP)))
                    currentpid += 1
                if specg is not None:
                    connection.execute('INSERT INTO photos(photo_id, event_id, photographer, description, date, file_extension) VALUES (?,?,?,?,?,?)',
                                       (currentpid, e1, 'Kate Allstadt', 'Spectrograms of velocity waveforms corrected for station sensitivity from closest stations with detections. Figure generation is automated and there may be artifacts present', date, os.path.join(relplacefigs, newfileSG)))
                    currentpid += 1
                if mt is not None:
                    connection.execute('INSERT INTO photos(photo_id, event_id, photographer, description, date, file_extension) VALUES (?,?,?,?,?,?)',
                                       (currentpid, e1, 'Kate Allstadt', 'Multitaper spectra (same units as Power Spectral Density) of velocity waveforms corrected for station sensitivity from closest stations with detections. Linear trend was removed and a taper was applied before generating spectra. Figure generation is automated and there may be artifacts present. Since station response was not removed to avoid adding artifacts at low frequencies for short and intermediate period stations, the lower frequency amplitudes should be interpreted with care.', date, os.path.join(relplacefigs, newfileSP)))
                    currentpid += 1
            print('updated database')
            connection.close()
        plt.close('all')


def grab_data(event_id, bufferperc=0.1, numstas=5, path='/Users/kallstadt/LSseis/landslideDatabase',
              detectHF=None, detectLP=None, both=False, minradius=0., maxradius=None, Zonly=False, reloadfile=False,
              loadfromfile=False, savedat=False, folderdat='data', merge=True, pad=True, fill_value=0., detrend='demean',
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
        if both:
            if v['detect_HF'] == 1 or v['detect_LP'] == 1:
                keep.append(v)
        elif detectHF is None and detectLP is None:
            keep.append(v)
        elif detectHF is True and detectLP is False:
            if v['detect_HF'] == 1 and v['detect_LP'] == 0:
                keep.append(v)
        elif detectHF is None and detectLP is True:
            if v['detect_LP'] == 1:
                keep.append(v)
        elif detectLP is True and detectHF is False:
            if v['detect_LP'] == 1 and v['detect_HF'] == 0:
                keep.append(v)
        elif detectLP is None and detectHF is True:
            if v['detect_HF'] == 1:
                keep.append(v)

    # Then keep only the numchans closest stations (all channels)
    if len(keep) == 0:
        return Stream()
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
        stsac = Stream()
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
                    stsac += reviewData.getdata_sac(newfilenames, attach_response=True, starttime=evDict['StartTime']-buffer_sec,
                                                    endtime=evDict['EndTime']+buffer_sec, chanuse=chanuse, savedat=False,
                                                    merge=merge, pad=pad, fill_value=fill_value)
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
                    stsac += reviewData.getdata_sac(filenames, attach_response=True, starttime=evDict['StartTime']-buffer_sec,
                                                    endtime=evDict['EndTime']+buffer_sec, chanuse=chanuse, savedat=False,
                                                    merge=merge, pad=pad, fill_value=fill_value)
                stsac = Stream([trace for trace in stsac if trace.max() != 0.0])  # Get rid of any empty ones
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
                                                   attach_response=True, clientname='IRIS', savedat=savedat,
                                                   folderdat=folderdat, filenamepref='iris_', reloadfile=reloadfile,
                                                   loadfromfile=loadfromfile, merge=merge, pad=pad, fill_value=fill_value)
            except Exception as e:
                print(e)
    if 'NCEDC' in evDict['DatLocation'] or 'NCEDC' in datlocs:
        stalist = [(k['Name'], k['Channel'], k['Network'], '*') for k in staDict2 if 'NCEDC' in k['source']]
        if len(stalist) != 0:
            sttemp += reviewData.getdata_exact(stalist, evDict['StartTime'] - buffer_sec, evDict['EndTime'] + buffer_sec,
                                               attach_response=True, clientname='NCEDC', savedat=savedat,
                                               folderdat=folderdat, filenamepref='ncedc_', reloadfile=reloadfile,
                                               loadfromfile=loadfromfile, merge=merge, pad=pad, fill_value=fill_value)

    # Delete any that are not in list
    nets, stas, chans = zip(*[[k['Network'], k['Name'], k['Channel']] for k in staDict])  # Exclude location code because inconsistencies here can cause data to not be found
    st = Stream()
    for n1, s1, c1 in zip(nets, stas, chans):
        st += sttemp.select(station=s1, network=n1, channel=c1)
    # Attach distaz
    st = findsta.attach_distaz(st, evDict['Latitude'], evDict['Longitude'], database=database)
    st = st.sort(keys=['rdist', 'channel'])

    # Preprocess
    if detrend is not None:
        st.detrend(detrend)

    return st


def make_figures(st, bufferperc=None, raw=True, HF=True, LP=True, timeseries=True, spectrograms=True, spectra=True,
                 HFlims=(1., 5.), HFoutput='VEL', LPlims=(20., 60.), LPoutput='DISP', taper=None, detrend='demean',
                 startline=None, endline=None, speccorr=None, removeoutliers=False, displacementspec=False):
    """
    Make static figures for each event
    :param HFlims: tuple or list of lower and upper frequency limits in Hz for HF filtering (1-5 Hz standard)
    :param HFoutput: Output type for HF station correction (VEL standard)
    :param taper: if None, no tapers, if 'bufferperc', will use bufferperc/2, if a float, will use that as the percent taper
    :param startline: UTCDateTime for start time of the event in the database, will add vertical line here in time series
    :param endline: UTCDateTime for end time of the event in the database, will add vertical line here in time series
    UPDATE THESE:
    spectra = multitaper, equiv to PSD
    """
    # build cosine filter that will be used
    cosfiltHF = (0.5*HFlims[0], HFlims[0], HFlims[1], 2*HFlims[1])
    cosfiltLP = (0.5/LPlims[1], 1/LPlims[1], 1/LPlims[0], 2/LPlims[0])

    dotaper = False
    if bufferperc is None and taper is not None:
        taper = taper
        dotaper = True
        if taper == 'bufferperc':
            taper = 0.05

    elif taper == 'bufferperc':
        if bufferperc is not None:
            taper = bufferperc/2.
        else:
            taper = 0.05
        dotaper = True

    if detrend is not None:
        st.detrend(detrend)

    removelater = Stream()

    # Do corrections for LP, HF
    if LP:
        stLP = st.copy()
        try:
            stLP.remove_response(output=LPoutput, pre_filt=cosfiltLP, taper=dotaper, taper_fraction=taper)
        except:
            temp = Stream()
            for i, trace in enumerate(stLP):
                try:
                    trace.remove_response(output=LPoutput, pre_filt=cosfiltLP)
                    temp = temp + trace
                except:
                    print 'Failed to remove response for %s, deleting this station' % (trace.stats.station + trace.stats.channel,)
                    removelater += trace
            stLP = temp.copy()
        if removeoutliers:
            for trace in stLP:
                if np.abs(trace.max()) > 100*np.median(np.abs(stLP.max())):
                    stLP.remove(trace)
                    removelater += trace

    if HF:
        stHF = st.copy()
        try:
            stHF.remove_response(output=HFoutput, pre_filt=cosfiltHF, taper=dotaper, taper_fraction=taper)
        except:
            temp = Stream()
            for i, trace in enumerate(stHF):
                try:
                    trace.remove_response(output=HFoutput, pre_filt=cosfiltHF)
                    temp = temp + trace
                except:
                    print 'Failed to remove response for %s, deleting this station' % (trace.stats.station + trace.stats.channel,)
                    removelater += trace
            stHF = temp.copy()
        if removeoutliers:
            for trace in stHF:
                print('made it here')
                if np.abs(trace.max()) > 10*np.median(np.abs(stHF.max())):
                    stHF.remove(trace)
                    removelater += trace

    # Figure out the colors, make sure same station has same color
    cycle = ['r', 'b', 'g', 'm', 'k', 'r', 'b', 'g', 'm', 'k', 'r', 'b', 'g', 'm', 'k', 'r', 'b', 'g', 'm', 'k', 'r', 'b', 'g', 'm', 'k', 'r', 'b', 'g', 'm', 'k', 'r', 'b', 'g', 'm', 'k']

    vlines = []
    if startline is not None:
        vlines.append(startline)
    if endline is not None:
        vlines.append(endline)

    # Now make figures
    if raw:
        colors1 = []
        stas = reviewData.unique_list([trace.stats.station for trace in st])
        for trace in st:
            ind = stas.index(trace.stats.station)
            colors1 += cycle[ind]
        figraw = reviewData.recsec(st, norm=True, maxtraces=15, quickdraw=True, figsize=(13, 14), colors=colors1,
                                   labelsize=14, addscale=False, unitlabel=None, convert=1., labelquickdraw=False,
                                   vlines=vlines, pad=True)
    else:
        figraw = None
    if HF:
        colors1 = []
        stas = reviewData.unique_list([trace.stats.station for trace in stHF])
        for tr in stHF.select(channel='BDF'):
            stHF.remove(tr)
        for trace in stHF:
            ind = stas.index(trace.stats.station)
            colors1 += cycle[ind]
        figHF = reviewData.recsec(stHF, norm=False, maxtraces=15, quickdraw=True, figsize=(13, 14), colors=colors1,
                                  labelsize=14, addscale=True, unitlabel='m/s', convert=1., labelquickdraw=False,
                                  vlines=vlines, pad=True)
    else:
        figHF = None
    if LP:
        colors1 = []
        stas = reviewData.unique_list([trace.stats.station for trace in stLP])
        for trace in stLP:
            ind = stas.index(trace.stats.station)
            colors1 += cycle[ind]
        figLP = reviewData.recsec(stLP, norm=False, maxtraces=15, quickdraw=True, figsize=(13, 14), colors=colors1,
                                  labelsize=14, addscale=True, unitlabel='m', convert=1., labelquickdraw=False,
                                  vlines=vlines, pad=True)
    else:
        figLP = None

    # Make spectrograms and spectra
    if spectrograms or spectra:
        stspec = st.copy()
        for trace in removelater:
            try:
                stspec.remove(trace)
            except:
                pass
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
        if removeoutliers:
            for trace in stspec:
                if np.abs(trace.max()) > 100*np.median(np.abs(stspec.max())):
                    stspec.remove(trace)

    if spectrograms:
        specg = reviewData.make_spectrogram(stspec, log1=True, maxtraces=15, maxPower=None, minPower=None, freqmax=None,
                                            labelsize=14, render=False)
    else:
        specg = None
    if spectra:
        if speccorr is None:
            colors1 = []
            stas = reviewData.unique_list([trace.stats.station for trace in stspec])
            for trace in stspec:
                ind = stas.index(trace.stats.station)
                colors1 += cycle[ind]
            j1, j2, mt = reviewData.make_multitaper(stspec, number_of_tapers=None, time_bandwidth=4., sine=False, recsec=True, colors1=colors1,
                                                    logx=True, logy=True, xunits='Hz', xlim=[0.001, 1000.], yunits='$m^2/s$', ylim=(10.**-22, 10.**-9),
                                                    render=False, detrend=detrend)  # ylim=(10.**-3, 10.**9)
        else:
            if speccorr == 'DISP':
                yunits = '$m/s$'
            elif speccorr == 'ACC':
                yunits = '$m^3/s$'
            else:
                yunits = '$m^2/s$'
            stspecD = st.copy()
            for trace in removelater:
                try:
                    stspecD.remove(trace)
                except:
                    pass
            try:
                stspecD.remove_response(output=speccorr, taper=dotaper, taper_fraction=taper)
            except:
                temp = Stream()
                for i, trace in enumerate(stspecD):
                    try:
                        trace.remove_response(output=speccorr, taper=dotaper, taper_fraction=taper)
                        temp = temp + trace
                    except:
                        print 'Failed to remove response for %s, deleting this station' % (trace.stats.station + trace.stats.channel,)
            if removeoutliers:
                for trace in stspecD:
                    if np.abs(trace.max()) > 100*np.median(np.abs(stspecD.max())):
                        stspecD.remove(trace)
            colors1 = []
            stas = reviewData.unique_list([trace.stats.station for trace in stspecD])
            for trace in stspecD:
                ind = stas.index(trace.stats.station)
                colors1 += cycle[ind]
            j1, j2, mt = reviewData.make_multitaper(stspecD, number_of_tapers=None, time_bandwidth=4., sine=False, recsec=True, colors1=colors1,
                                                    logx=True, logy=True, xunits='Hz', xlim=[0.001, 1000.], yunits=yunits,
                                                    render=False, ylim=(10.**-25, 10.**-9))
    else:
        mt = None

    return figraw, figHF, figLP, specg, mt
