#!/usr/bin/env python

"""
Functions for populating and updating the database
"""

import findsta
from reviewData import reviewData
import numpy as np
import sqlite3 as lite
from obspy import Stream, UTCDateTime
from obspy.clients.fdsn import Client as FDSN_Client
import glob
import urllib.request
import os
from sigproc import sigproc
import shutil
import datetime


def initial_populate(event_ids, minradius=0., maxradius=500., clients=['IRIS'],
                     database='/Users/kallstadt/LSseis/landslideDatabase/lsseis.db'):
    """Performs initial population of stations and sta_nearby tables for a single event id or list of event ids.

    Args:
        event_ids: either single integer or list or numpy array of integers, eg np.arange(62,65)
        maxradius (float): Maximum radius to go out to in searching for nearby stations
        clients (list): list of FDSN clients to search, see options here: https://docs.obspy.org/packages/obspy.clients.fdsn.html
        database (str): file path to sqlite3 database
    """
    if type(event_ids) is int:
        event_ids = [event_ids]
    for event_id in event_ids:
        event_id = int(event_id)
        evdict = findsta.getEventInfo(event_id, database=database)
        if 'AlaskaData' in evdict['DatLocation']:  # Need to trick the system because IRIS doesn't have the actual install dates for these stations
            stb4 = UTCDateTime(2009, 9, 1, 0, 0, 0)
        else:
            stb4 = evdict['StartTime']
        for client in clients:

            inventory = reviewData.get_stations(evdict['Latitude'], evdict['Longitude'], stb4,
                                                clients=[client], minradiuskm=minradius,
                                                maxradiuskm=maxradius, chan='BH?,EH?,HH?,EL?,HN?,EN?,CH?,DH?,LH?,SL?')

            #print(len(inventory.get_contents()['channels']))
            populate_station_tables(inventory, client=client, database=database)
            populate_station_event_table(event_id, inventory, database=database)


def populate_station_tables(inventory, client=None, database='/Users/kallstadt/LSseis/landslideDatabase/lsseis.db'):
    """Helper code used to put the station info into the stations table of database if it's not already there

    Args:
        inventory: Obspy Inventory object containing stations to try to add to table.
        client (str): FDSN client from which station information was obtained. Used to define source of station info
            in table.
        database (str): file path to sqlite3 database
    """
    if client is None:
        source = 'unknown'
    else:
        source = client
    #loop over each line, parse info, see if it already is in table, if not, insert it
    connection = None
    connection = lite.connect(database)
    with connection:
        cursor = connection.cursor()
        val = 0
        val2 = 0
        networks = inventory.networks
        for net in networks:
            stations = net.stations
            for sta in stations:
                channels = sta.channels
                for chan in channels:
                    #try:
                    #see if the station is already there
                    cursor_output = (cursor.execute(
                        'SELECT * FROM stations WHERE Name=? AND Channel=? AND Network=?', (sta.code, chan.code, net.code)))
                    #except Exception as e:
                    #    print(e)
                    retrieved_data = cursor_output.fetchall()
                    if len(retrieved_data) == 0:  # if the station isn't already there, put it there
                        cursor.execute(
                            'INSERT INTO stations(Network,Name,LocationCode,Channel,Latitude,Longitude,Elevation_masl,source) VALUES(?,?,?,?,?,?,?,?)',
                            (net.code, sta.code, chan.location_code, chan.code, chan.latitude, chan.longitude, chan.elevation, source))
                        val += 1
                    else:
                        val2 += 1
        print(('added %s entries to stations table, %s already were there' % (val, val2)))


def populate_station_event_table(event_id, inventory, database='/Users/kallstadt/LSseis/landslideDatabase/lsseis.db',
                                 update=True):
    """Helper code to put the station in the sta_nearby table for this event and calculate distance, az, baz

    Args:
        event_id (int): event id of landslide event
        inventory: Obspy Inventory object containing stations to try to add to table.
        database (str): file path to sqlite3 database
        update (bool): if True, will recalculate distances and update all station-event pairs that were already
            in the database.

    """
    #get event lat lon from database
    connection = None
    connection = lite.connect(database)
    with connection:
        cursor = connection.cursor()
        #try:
        cursor_output = cursor.execute('SELECT Latitude, Longitude FROM events WHERE Eid=?', (event_id,))
        #except Exception as e:
        #    print(e)
        #    return
    retrieved_data = cursor_output.fetchall()
    event_lat, event_lon = retrieved_data[0]

    #get Sid and lat lon for each entry in line
    val = 0
    val1 = 0
    val2 = 0
    valbad = 0

    networks = inventory.networks
    with connection:
        for net in networks:
            stations = net.stations
            for sta in stations:
                channels = sta.channels
                for chan in channels:
                    #Get station id from station table
                    #print('%s.%s.%s' % (sta.code, chan.code, net.code))
                    cursor_output = cursor.execute(
                        'SELECT Sid, Latitude, Longitude FROM stations WHERE Name=? AND Channel=? AND Network=?', (sta.code, chan.code, net.code))
                    retrieved_data = cursor_output.fetchone()
                    Sid, sta_lat, sta_lon = retrieved_data

                    # See if entry is already there
                    #try:
                    cursor_output = cursor.execute('SELECT SRid FROM sta_nearby WHERE event_id=? AND station_id=?',
                                                   (event_id, Sid))
                    retrieved_data = cursor_output.fetchall()
                    #except:
                    #    retrieved_data = []

                    if len(retrieved_data) == 0:
                        #if it isn't already there, calculate stuff and save it in the table
                        backazimuth, azimuth, distance = reviewData.pyproj_distaz(sta_lat, sta_lon, event_lat, event_lon)
                        #try:
                        cursor.execute('INSERT INTO sta_nearby(event_id,station_id,stasource_radius_km,az,baz) VALUES(?,?,?,?,?)',
                                       (event_id, Sid, '%3.3f' % distance, '%3.2f' % azimuth, '%3.2f' % backazimuth))
                        val += 1
                        #except Exception as f:
                        #    print(f)
                        #    valbad += 1
                    elif update is True:
                        #try:
                        #result = client.distaz(sta_lat, sta_lon, event_lat, event_lon)
                        backazimuth, azimuth, distance = reviewData.pyproj_distaz(sta_lat, sta_lon, event_lat, event_lon)
                        cursor.execute("UPDATE sta_nearby SET stasource_radius_km = ?, az = ?, baz = ? WHERE event_id =? AND station_id =?;",
                                       ('%3.3f' % distance, '%3.2f' % azimuth, '%3.2f' % backazimuth, event_id, Sid))
                        val1 += 1
                        #except Exception as f:
                        #    print(f)
                        #    valbad += 1
                    else:
                        val2 += 1
    print(('added %s entries to sta_nearby table, updated %s entries, %s left as is, %s not added because of error' % (val, val1, val2, valbad)))


def remove_events(event_ids, savecopy=None, gisfiles=False, photos=False, information=False,
                  database='/Users/kallstadt/LSseis/landslideDatabase/lsseis.db'):
    """Removes event_ids from event table and also removes any entries in sta_nearby table for this event.

    Args:
        event_ids: either single integer or list or numpy array of integers, eg np.arange(62,65)
        savecopy (str): full file name for a copy to save of database in case something goes wrong with .db extension
        gisfiles (bool): if True, remove corresponding entries in gisfiles table
        photos (bool): if True, remove corresponding entries in photos table
        information (bool): if True, remove corresponding entries in information table
        database (str): file path to sqlite3 database
    """
    if type(event_ids) is int:
        event_ids = [event_ids]
    if savecopy is not None:
        dir1 = savecopy.split('/')
        try:
            os.makedirs(os.path.join(*dir1[:-1]))
        except:
            pass
        try:
            time1 = datetime.datetime.utcnow().strftime('%d%b%Y_%H%M')
            filename = savecopy.split('.')[0] + '_' + time1 + '.' + savecopy.split('.')[-1]
            # add datetime to end
            shutil.copy(database, filename)
        except Exception as e:
            print(e)
            print('Unable to save a copy of the database, exiting')
            return

    connection = None
    connection = lite.connect(database)
    for eid in event_ids:
        cursor = connection.cursor()
        try:
            # delete the event table entry
            cursor.execute('DELETE FROM events WHERE Eid=?', (eid,))
            # delete sta_nearby entries
            cursor.execute('DELETE FROM sta_nearby WHERE event_id=?', (eid,))
            if gisfiles:
                cursor.execute('DELETE FROM gisfiles WHERE event_id=?', (eid,))
            if photos:
                cursor.execute('DELETE FROM photos WHERE event_id=?', (eid,))
            if information:
                cursor.execute('DELETE FROM information WHERE event_id=?', (eid,))
            connection.commit()
        except Exception as e:
            print(e)
            # roll back the changes if anything failed
            connection.rollback()


def recalculate_distances(event_ids, database='/Users/kallstadt/LSseis/landslideDatabase/lsseis.db'):
    """Recalculate source to station distances. Should be done if lat/lons were updated for either source or stations in the table

    Args:
        event_ids: either single integer or list or numpy array of integers, eg np.arange(62,65)
        database (str): file path to sqlite3 database
    """
    if type(event_ids) is int:
        event_ids = [event_ids]
    connection = None
    connection = lite.connect(database)

    for event_id in event_ids:
    # Get all the SRids and corresponding station_id's for event_id
        evdict = findsta.getEventInfo(event_id, database=database)
        cursor = connection.cursor()
        cursor_output = (cursor.execute('SELECT SRid, station_id FROM sta_nearby WHERE event_id=?', (event_id,)))
        retrieved_data = cursor_output.fetchall()
        SRid = [dat[0] for dat in retrieved_data]
        station_id = [dat[1] for dat in retrieved_data]
        for i, staid in enumerate(station_id):
            cursor_output = (cursor.execute('SELECT Latitude, Longitude FROM stations WHERE Sid=?', (staid,)))
            retrieved_data = cursor_output.fetchall()
            sta_lat = retrieved_data[0][0]
            sta_lon = retrieved_data[0][1]
            backazimuth, azimuth, distance = reviewData.pyproj_distaz(sta_lat, sta_lon, evdict['Latitude'],
                                                                      evdict['Longitude'])
            with connection:
                connection.execute("UPDATE sta_nearby SET stasource_radius_km = ?, az = ?, baz = ? WHERE SRid=?",
                                   ('%3.3f' % distance, '%3.2f' % azimuth, '%3.2f' % backazimuth, SRid[i]))
    connection.close()


def review_event(event_id, buffer_sec=100., minradius=0., maxradius=200., intincrkm=100.,
                 maxreachedHF=False, maxreachedLP=False, HFlims=(1., 5.), LPlims=(20., 60.),
                 LPoutput='DISP', maxtraces=15, taper=0.05, clients=['IRIS'],
                 database='/Users/kallstadt/LSseis/landslideDatabase/lsseis.db',
                 path='/Users/kallstadt/LSseis/landslideDatabase'):
    """
    Initial review of seismic data for a given event, pulls data from IRIS and other sources and
    displays in HF and LP bands (depending on input choices), user should then delete all traces
    in each where the signal is not visible to the human eye, leaving only the ones that should
    be considered "detections" behind. Then quit interactive plotting by pressing q and follow
    prompts. Will populate database out to maxradius and maxradius + n*intincrkm.
    Note: This must be done in ipython with pylab enabled for interactive plotting to work. May also require
    that an interactive backend be used for matplotlib (Qt5Agg works).

    Args:
        event_id (int): id of event to review
        buffer_sec (float): Number of seconds to add on either end of start and end time (makes viewing easier)
        minradius (float): minimum distance in km to search for stations (from source)
        maxradius (float): maximum distance in km to search for stations (from source)
        intincrkm (float): increment to increase by when expanding search outward past maxradius
        maxreachedHF (bool): Has the maximum extent of the high frequencies already been reached in a previous
            run of review data and saved in database? (will skip analysis of high frequencies)
        maxreachedLP (bool): Has the maximum extent of the long periods already been reached in a previous
            run of review data and saved in database? (will skip analysis of high frequencies)
        HFlims: tuple or list of lower and upper frequency limits in Hz for HF filtering (1-5 Hz standard)
        LPlims: tuple or list of lower and upper period limits in seconds for LP station correction (in DISP),
            (20., 60.) standard
        LPoutput (str): Output type for LP station correction ('DISP' standard, could also be 'VEL' or 'ACC')
        maxtraces (int): Number of traces to view at a time
        taper (0.05): percentage taper (or None for no tapering)
        database (str): Full file path of database file location
        path (str): Path to location of sac files (upstream from relative file paths listed in database)

    """
    evDict = findsta.getEventInfo(event_id, database=database)
    print((('Now analysing Eid %s - %s') % (event_id, evDict['Name'])))
    # Populate database with stations, if it isn't already populated, out to maxradius
    staDict = findsta.getStaInfo(event_id, database=database)
    dists = [sdct['stasource_radius_km'] for sdct in staDict]
    if np.max(dists) < maxradius:
        print ('database not fully populated to maxradius, populating now')
        initial_populate(event_id, minradius=np.max(dists), maxradius=maxradius, clients=clients,
                         database=database)

    # Load data from station list to maxradius
    if maxreachedHF is True and maxreachedLP is False:  # If only looking at LP, don't bother downloading other stuff
        staDict = findsta.getStaInfo(event_id, maxradius=maxradius, minradius=minradius,
                                     chanuse='BHZ,BHE,BHN,BH1,BH2,HHZ,HHE,HHN,HH1,HH2,LHZ,LHN,LHE,LH1,LH2,LH3,SLZ,SLN,SLE', database=database)
    else:
        staDict = findsta.getStaInfo(event_id, maxradius=maxradius, minradius=minradius, database=database)
    datlocs = reviewData.unique_list([sdct['source'] for sdct in staDict])

    if len(staDict) > 0:
        # Download data from their respective sources
        st = Stream()
        if 'IRIS' in evDict['DatLocation'] or 'IRIS' in datlocs:
            try:
                stalist, netlist, chanlist = list(zip(*[[sdct['Name'], sdct['Network'],
                                                  sdct['Channel']] for sdct in staDict if 'IRIS' in sdct['source']]))
                st += reviewData.getdata(','.join(reviewData.unique_list(netlist)), ','.join(reviewData.unique_list(stalist)),
                                         '*', ','.join(reviewData.unique_list(chanlist)), evDict['StartTime']-buffer_sec,
                                         evDict['EndTime']+buffer_sec, savedat=False, detrend='demean')
            except Exception as e:
                print(e)
        if 'NCEDC' in evDict['DatLocation'] or 'NCEDC' in datlocs:
            try:
                stalist, netlist, chanlist = list(zip(*[[sdct['Name'], sdct['Network'],
                                                  sdct['Channel']] for sdct in staDict if 'NCEDC' in sdct['source']]))
                st += reviewData.getdata(','.join(reviewData.unique_list(netlist)), ','.join(reviewData.unique_list(stalist)),
                                         '*', ','.join(reviewData.unique_list(chanlist)), evDict['StartTime']-buffer_sec,
                                         evDict['EndTime']+buffer_sec, savedat=False, clientname='NCEDC', detrend='demean')
            except Exception as e:
                print(e)
        if evDict['DatLocation'] is None:
            print('You need to populate the DatLocation field for this event, no sac files loaded')
        if 'sac' in evDict['DatLocation']:
            datloc1 = evDict['DatLocation'].split(',')
            datloc1 = [x.strip() for x in datloc1 if 'sac' in x]
            for datl in datloc1:
                fullpath = os.path.join(path, datl.split(':')[1])
                filenames = glob.glob(fullpath)
                if len(filenames) > 0:
                    if 'AlaskaData' in datl:  # Need to do some cheating to attach response info if Iliamna sac data - attach oldest response info available at IRIS for each station
                        # Only keep filenames of stations we want to read in
                        newfilenames = []
                        namelist = [sdct['Name'] for sdct in staDict]
                        for filen in filenames:
                            nam = filen.split('/')[-1].split('.')[0]
                            if nam in namelist:
                                newfilenames.append(filen)
                        stsac = reviewData.getdata_sac(newfilenames, attach_response=True, starttime=evDict['StartTime']-buffer_sec,
                                                       endtime=evDict['EndTime']+buffer_sec, detrend='demean')
                        stsac = Stream([trace for trace in stsac if trace.max() != 0.0])  # Get rid of any empty ones
                        originalstt = []
                        # Get earliest start time from IRIS for AK and AV
                        url1 = ('http://service.iris.edu/fdsnws/station/1/query?network=AV&level=station&format=text&nodata=404')
                        url2 = ('http://service.iris.edu/fdsnws/station/1/query?network=AK&level=station&format=text&nodata=404')
                        f = urllib.request.urlopen(url1)
                        file1 = f.read()
                        lines1 = [line.split('|') for line in file1.split('\n')[1:]]
                        stdict = {}
                        f.close()
                        f = urllib.request.urlopen(url2)
                        file2 = f.read()
                        lines2 = [line.split('|') for line in file2.split('\n')[1:]]
                        f.close()
                        for line in lines1[:-1]:
                            stdict[line[1]] = (line[0], UTCDateTime(line[6]))
                        for line in lines2[:-1]:
                            if line[1] in list(stdict.keys()):
                                # Take the minimum if shows up for both network codes
                                stdict[line[1]] = (line[0], np.min([stdict[line[1]], UTCDateTime(line[6])]))
                            else:
                                stdict[line[1]] = (line[0], UTCDateTime(line[6]))
                        for temp in stsac:
                            originalstt.append(temp.stats.starttime)
                            if temp.stats.channel[0] == 'S':
                                temp.stats.channel = 'E'+temp.stats.channel[1:]
                            # Make sure network code is consistent with IRIS and delete events with nothing returned from IRIS
                            if temp.stats.station in stdict:
                                temp.stats.network = stdict[temp.stats.station][0]
                            else:
                                print('could not attach response info for %s, station correction will not work' % temp.stats.station)
                                continue
                            if 'response' not in list(temp.stats.keys()):
                                # Change starttime of temp to earliest start time at IRIS for this station
                                try:
                                    temp.stats.starttime = stdict[temp.stats.station][1] + 86400.
                                except:
                                    print('could not attach response info for %s, station correction will not work' % temp.stats.station)
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
                                                       endtime=evDict['EndTime']+buffer_sec, detrend='demean')
                    st += stsac

        # Add distances
        st = findsta.attach_distaz(st, evDict['Latitude'], evDict['Longitude'], database=database)
        st = st.sort(keys=['rdist', 'channel'])
        # Get rid of sac files etc. that are outside current range
        st = Stream([trace for trace in st if trace.stats.rdist < maxradius and trace.stats.rdist > minradius])
        if taper is not None:
            st.taper(max_percentage=taper, type='cosine')

        if maxreachedHF is False:
            # Split into hf and lp streams and preprocess
            st_hf = st.copy()
            # pre-filter to hf_range
            st_hf.filter('bandpass', freqmin=HFlims[0], freqmax=HFlims[1])
            # pre-correct to lp range, delete E* channels and any that won't station correct
            print ('Interactive review of high frequencies starting\n')
            # Interactive review of high frequencies
            zp = reviewData.InteractivePlot(st_hf, maxtraces=maxtraces, textline=['HF filtered data for your review',
                                            'Delete bad traces and take note of distance when you reach maximum observation',
                                            'hit Q to exit when ready to continue'])
            # Get info about what was kept and what wasn't
            st_hfgood = zp.st_current
            hfbad = zp.deleted

            connection = None
            connection = lite.connect(database)

            # See if maxradiusHF was reached
            maxdist = np.max([trace.stats.rdist for trace in st_hfgood])
            feedback = input(('is %4.1f km the maximum distance observed for HF? Y or N\n') % (maxdist,))
            if feedback.lower() != 'n' and feedback.lower() != 'y':
                feedback = input(('Try again, is %4.1f km the maximum distance observed for HF? Y or N\n') % (maxdist,))
            if feedback.lower() == 'n':
                feedback = input('was the maximum distance less than %4.1f km?\n' % (maxdist,))
                if feedback.lower() == 'y':
                    maxdist = float(input('Enter the maximum distance observed in km (round up to nearest km)\n'))
                    # delete all stations that were greater than maximum distance and that are in st_hfbad, and put good stations in database
                    goodstas = [str(trace.stats.station) for trace in st_hfgood if trace.stats.rdist <= maxdist]
                    goodchans = [str(trace.stats.channel) for trace in st_hfgood if trace.stats.rdist <= maxdist]
                    goodnets = [str(trace.stats.network) for trace in st_hfgood if trace.stats.rdist <= maxdist]
                    badstas = [str(trace.stats.station) for trace in st_hfgood if trace.stats.rdist > maxdist] + [str(g.split('.')[0]) for g in hfbad]
                    badchans = [str(trace.stats.channel) for trace in st_hfgood if trace.stats.rdist > maxdist] + [str(g.split('.')[1]) for g in hfbad]
                    badnets = [str(trace.stats.network) for trace in st_hfgood if trace.stats.rdist > maxdist] + [str(g.split('.')[3].split(' - ')[0]) for g in hfbad]
                    maxreachedHF = True
                    # Insert maxdisthf into database
                    with connection:
                        cursor = connection.cursor()
                        try:
                            cursor.execute('UPDATE events SET maxdistHF_km = ? WHERE Eid = ?', (maxdist, event_id))
                        except Exception as e:
                            print(e)
                if feedback.lower() == 'n':
                    print('Filling current info into database, then will load more distant data for further analysis\n')
                    badstas = [str(g.split('.')[0]) for g in hfbad]
                    badchans = [str(g.split('.')[1]) for g in hfbad]
                    badnets = [str(g.split('.')[3].split(' - ')[0]) for g in hfbad]
                    goodstas = [str(trace.stats.station) for trace in st_hfgood]
                    goodchans = [str(trace.stats.channel) for trace in st_hfgood]
                    goodnets = [str(trace.stats.network) for trace in st_hfgood]
            elif feedback.lower() == 'y':
                # put info in database
                badstas = [str(g.split('.')[0]) for g in hfbad]
                badchans = [str(g.split('.')[1]) for g in hfbad]
                badnets = [str(g.split('.')[3].split(' - ')[0]) for g in hfbad]
                goodstas = [str(trace.stats.station) for trace in st_hfgood]
                goodchans = [str(trace.stats.channel) for trace in st_hfgood]
                goodnets = [str(trace.stats.network) for trace in st_hfgood]
                maxreachedHF = True
                # Insert maxdisthf into database
                with connection:
                    cursor = connection.cursor()
                    try:
                        cursor.execute('UPDATE events SET maxdistHF_km = ? WHERE Eid = ?', (maxdist, event_id))
                    except Exception as e:
                        print(e)
            # Fill info into database
            Sids_bad = []
            Sids_good = []
            with connection:
                cursor = connection.cursor()
                for i, badsta in enumerate(badstas):
                    try:
                        cursor_output = cursor.execute('SELECT Sid FROM stations WHERE Name=? AND Network=? AND Channel=?',
                                                       (badstas[i], badnets[i], badchans[i]))
                        retrieved_data = cursor_output.fetchall()
                        if retrieved_data:
                            Sids_bad += [val[0] for val in retrieved_data]
                    except Exception as e:
                        print(e)
                for i, goodsta in enumerate(goodstas):
                    try:
                        cursor_output = cursor.execute('SELECT Sid FROM stations WHERE Name=? AND Network=? AND Channel=?',
                                                       (goodstas[i], goodnets[i], goodchans[i]))
                        retrieved_data = cursor_output.fetchall()
                        if retrieved_data:
                            Sids_good += [val[0] for val in retrieved_data]
                    except Exception as e:
                        print(e)
                if Sids_bad:
                    for Sid in Sids_bad:
                        cursor.execute('UPDATE sta_nearby SET detect_HF = ? WHERE station_id = ? AND event_id = ?',
                                       (0, Sid, event_id))
                if Sids_good:
                    for Sid in Sids_good:
                        cursor.execute('UPDATE sta_nearby SET detect_HF = ? WHERE station_id = ? AND event_id = ?',
                                       (1, Sid, event_id))
                # if maxreachedHF is True:
                #     cursor.execute('UPDATE sta_nearby SET detect_HF = ? WHERE stasource_radius_km > ? AND event_id = ?',
                #                    (0, float(maxdist), event_id))

        if maxreachedLP is False:
            st_lp = st.copy()
            st_lp = st_lp.select(channel='B*') + st_lp.select(channel='H*') + st_lp.select(channel='L*')
            if len(st_lp) > 0:
                #reorder by distance
                st_lp = st_lp.sort(keys=['rdist', 'channel'])
                temp = Stream()
                cosfilt = (0.5/LPlims[1], 1/LPlims[1], 1/LPlims[0], 2/LPlims[0])
                for i, trace in enumerate(st_lp):
                    try:
                        trace.remove_response(output=LPoutput, pre_filt=cosfilt)
                        temp = temp + trace
                    except:
                        print('Failed to remove response for %s, deleting this station' % (trace.stats.station + trace.stats.channel,))
                st_lp = temp

                print ('Interactive review of long periods starting\n')
                # Interactive review of long periods
                zp = reviewData.InteractivePlot(st_lp, maxtraces=maxtraces, textline=['Station corrected LP data for your review',
                                                'Delete bad traces and take note of distance when you reach maximum observation',
                                                'hit Q to exit when ready to continue'])
                # Get info about what was kept and what wasn't
                st_lpgood = zp.st_current
                lpbad = zp.deleted

                connection = None
                connection = lite.connect(database)

                # See if maxradiuslp was reached
                maxdist = np.max([trace.stats.rdist for trace in st_lpgood])
                feedback = input(('is %4.1f km the maximum distance observed for lp? Y or N\n') % (maxdist,))
                if feedback.lower() == 'n':
                    feedback = input('was the maximum distance less than %4.1f km?\n' % (maxdist,))
                    if feedback.lower() == 'y':
                        maxdist = float(input('Enter the maximum distance observed in km (round up to nearest km)\n'))
                        # delete all stations that were greater than maximum distance and that are in st_lpbad, and put good stations in database
                        goodstas = [str(trace.stats.station) for trace in st_lpgood if trace.stats.rdist < maxdist]
                        goodchans = [str(trace.stats.channel) for trace in st_lpgood if trace.stats.rdist < maxdist]
                        goodnets = [str(trace.stats.network) for trace in st_lpgood if trace.stats.rdist < maxdist]
                        badstas = [str(trace.stats.station) for trace in st_lpgood if trace.stats.rdist > maxdist] + [str(g.split('.')[0]) for g in lpbad]
                        badchans = [str(trace.stats.channel) for trace in st_lpgood if trace.stats.rdist > maxdist] + [str(g.split('.')[1]) for g in lpbad]
                        badnets = [str(trace.stats.network) for trace in st_lpgood if trace.stats.rdist > maxdist] + [str(g.split('.')[3].split(' - ')[0]) for g in lpbad]
                        maxreachedLP = True
                        # Insert maxdistlp into database
                        with connection:
                            cursor = connection.cursor()
                            try:
                                cursor.execute('UPDATE events SET maxdistLP_km = ? WHERE Eid = ?', (maxdist, event_id))
                            except Exception as e:
                                print(e)
                    if feedback.lower() == 'n':
                        print('Filling current info into database, then will load more distant data for further analysis\n')
                        badstas = [str(g.split('.')[0]) for g in lpbad]
                        badchans = [str(g.split('.')[1]) for g in lpbad]
                        badnets = [str(g.split('.')[3].split(' - ')[0]) for g in lpbad]
                        goodstas = [str(trace.stats.station) for trace in st_lpgood]
                        goodchans = [str(trace.stats.channel) for trace in st_lpgood]
                        goodnets = [str(trace.stats.network) for trace in st_lpgood]
                elif feedback.lower() == 'y':
                    # put info in database
                    badstas = [str(g.split('.')[0]) for g in lpbad]
                    badchans = [str(g.split('.')[1]) for g in lpbad]
                    badnets = [str(g.split('.')[3].split(' - ')[0]) for g in lpbad]
                    goodstas = [str(trace.stats.station) for trace in st_lpgood]
                    goodchans = [str(trace.stats.channel) for trace in st_lpgood]
                    goodnets = [str(trace.stats.network) for trace in st_lpgood]
                    maxreachedLP = True
                    # Insert maxdistlp into database
                    with connection:
                        cursor = connection.cursor()
                        try:
                            cursor.execute('UPDATE events SET maxdistLP_km = ? WHERE Eid = ?', (maxdist, event_id))
                        except Exception as e:
                            print(e)
                # Fill info into database
                Sids_bad = []
                Sids_good = []
                with connection:
                    cursor = connection.cursor()
                    for i, badsta in enumerate(badstas):
                        try:
                            cursor_output = cursor.execute('SELECT Sid FROM stations WHERE Name=? AND Network=? AND Channel=?',
                                                           (badstas[i], badnets[i], badchans[i]))
                            retrieved_data = cursor_output.fetchall()
                            if retrieved_data:
                                Sids_bad += [val[0] for val in retrieved_data]
                        except Exception as e:
                            print(e)
                    for i, goodsta in enumerate(goodstas):
                        try:
                            cursor_output = cursor.execute('SELECT Sid FROM stations WHERE Name=? AND Network=? AND Channel=?',
                                                           (goodstas[i], goodnets[i], goodchans[i]))
                            retrieved_data = cursor_output.fetchall()
                            if retrieved_data:
                                Sids_good += [val[0] for val in retrieved_data]
                        except Exception as e:
                            print(e)
                    if Sids_bad:
                        for Sid in Sids_bad:
                            cursor.execute('UPDATE sta_nearby SET detect_LP = ? WHERE station_id = ? AND event_id = ?',
                                           (0, Sid, event_id))
                    if Sids_good:
                        for Sid in Sids_good:
                            cursor.execute('UPDATE sta_nearby SET detect_LP = ? WHERE station_id = ? AND event_id = ?',
                                           (1, Sid, event_id))
                    # if maxreachedLP is True:
                    #     cursor.execute('UPDATE sta_nearby SET detect_LP = ? WHERE stasource_radius_km > ? AND event_id = ?',
                    #                    (0, float(maxdist), event_id))
            else:
                print('no stations in current range')

    # Continue searching for maximum distance
    while maxreachedHF is False or maxreachedLP is False:
        print('Going out another %3.0f km' % (intincrkm,))
        newmin = maxradius
        maxradius = maxradius + intincrkm

        # Populate database with stations, if it isn't already populated, out to maxradius
        staDict = findsta.getStaInfo(event_id, database=database)
        dists = [sdct['stasource_radius_km'] for sdct in staDict]
        if np.max(dists) < maxradius:
            print ('database not fully populated to maxradius, populating now')
            initial_populate(event_id, minradius=newmin, maxradius=maxradius, clients=clients,
                             database=database)

        # Load data from station list to maxradius
        if maxreachedHF is True and maxreachedLP is False:  # If only looking at LP, don't bother downloading other stuff
            staDict = findsta.getStaInfo(event_id, maxradius=maxradius, minradius=newmin,
                                         chanuse='BHZ,BHE,BHN,BH1,BH2,HHZ,HHE,HHN,HH1,HH2,LHZ,LHN,LHE,LH1,LH2,LH3,SLZ,SLN,SLE', database=database)
        else:
            staDict = findsta.getStaInfo(event_id, maxradius=maxradius, minradius=newmin, database=database)
        datlocs = reviewData.unique_list([sdct['source'] for sdct in staDict])

        if len(staDict) == 0:
            print ('no stations in given radius')
            continue

        # Load and preprocess data - TO DO: ONLY LOAD STATIONS THAT HAVENT BEEN REVIEWED YET?
        st = Stream()

        if 'IRIS' in evDict['DatLocation'] or 'IRIS' in datlocs:
            try:
                stalist, netlist, chanlist = list(zip(*[[sdct['Name'], sdct['Network'], sdct['Channel']] for sdct in staDict if 'IRIS' in sdct['source']]))
                # Remove any BDF's
                #stalist, netlist, chanlist = zip(*[[stalist[k], netlist[k], chanlist[k]] for k in range(len(chanlist)) if 'BDF' not in chanlist[k]])
                st += reviewData.getdata(','.join(reviewData.unique_list(netlist)), ','.join(reviewData.unique_list(stalist)),
                                         '*', ','.join(reviewData.unique_list(chanlist)), evDict['StartTime']-buffer_sec,
                                         evDict['EndTime']+buffer_sec, savedat=False)
            except:
                print('No data from IRIS found in current distance interval')
        if 'NCEDC' in evDict['DatLocation'] or 'NCEDC' in datlocs:
            try:
                stalist, netlist, chanlist = list(zip(*[[sdct['Name'], sdct['Network'], sdct['Channel']] for sdct in staDict if 'NCEDC' in sdct['source']]))
                st += reviewData.getdata(','.join(reviewData.unique_list(netlist)), ','.join(reviewData.unique_list(stalist)),
                                         '*', ','.join(reviewData.unique_list(chanlist)), evDict['StartTime']-buffer_sec,
                                         evDict['EndTime']+buffer_sec, savedat=False, clientname='NCEDC')
            except:
                print('No data from NCEDC found in current distance interval')
        if evDict['DatLocation'] is None:
            print('You need to populate the DatLocation field for this event, no sac files loaded')
        elif 'sac' in evDict['DatLocation']:
            st += stsac  # already loaded in above
        if len(st) == 0:
            continue

        # Add distances
        st = findsta.attach_distaz(st, evDict['Latitude'], evDict['Longitude'], database=database)
        st = st.sort(keys=['rdist', 'channel'])
        # Get rid of sac files etc. that are outside current range
        st = Stream([trace for trace in st if trace.stats.rdist < maxradius and trace.stats.rdist > newmin])
        if taper is not None:
            st.taper(max_percentage=taper, type='cosine')

        if maxreachedHF is False:
            st_hf = st.copy()
            # pre-filter to hf_range
            st_hf.filter('bandpass', freqmin=HFlims[0], freqmax=HFlims[1])
            # pre-correct to lp range, delete E* channels and any that won't station correct
            print ('Interactive review of high frequencies starting\n')
            # Interactive review of high frequencies
            zp = reviewData.InteractivePlot(st_hf, maxtraces=maxtraces, textline=['HF filtered data for your review',
                                            'Delete bad traces and take note of distance when you reach maximum observation',
                                            'hit Q to exit when ready to continue'])
            #zp.connect()
            # Get info about what was kept and what wasn't
            st_hfgood = zp.st_current
            hfbad = zp.deleted

            connection = None
            connection = lite.connect(database)
            maxdist = np.max([trace.stats.rdist for trace in st_hfgood])
            feedback = input(('is %4.1f km the maximum distance observed for HF? Y or N\n') % (maxdist,))
            if feedback.lower() != 'n' and feedback.lower != 'y':
                feedback = input(('Try again, is %4.1f km the maximum distance observed for HF? Y or N\n') % (maxdist,))
            if feedback.lower() == 'n':
                feedback = input('was the maximum distance less than %4.1f km?\n' % (maxdist,))
                if feedback.lower() == 'y':
                    maxdist = float(input('Enter the maximum distance observed in km (round up to nearest km)\n'))
                    # delete all stations that were greater than maximum distance and that are in st_hfbad, and put good stations in database
                    goodstas = [str(trace.stats.station) for trace in st_hfgood if trace.stats.rdist <= maxdist]
                    goodchans = [str(trace.stats.channel) for trace in st_hfgood if trace.stats.rdist <= maxdist]
                    goodnets = [str(trace.stats.network) for trace in st_hfgood if trace.stats.rdist <= maxdist]
                    badstas = [str(trace.stats.station) for trace in st_hfgood if trace.stats.rdist > maxdist] + [str(g.split('.')[0]) for g in hfbad]
                    badchans = [str(trace.stats.channel) for trace in st_hfgood if trace.stats.rdist > maxdist] + [str(g.split('.')[1]) for g in hfbad]
                    badnets = [str(trace.stats.network) for trace in st_hfgood if trace.stats.rdist > maxdist] + [str(g.split('.')[3].split(' - ')[0]) for g in hfbad]
                    maxreachedHF = True
                    # Insert maxdisthf into database
                    with connection:
                        cursor = connection.cursor()
                        try:
                            cursor.execute('UPDATE events SET maxdistHF_km = ? WHERE Eid = ?', (maxdist, event_id))
                        except Exception as e:
                            print(e)
                if feedback.lower() == 'n':
                    print('Filling current info into database, then will load more distant data for further analysis\n')
                    badstas = [str(g.split('.')[0]) for g in hfbad]
                    badchans = [str(g.split('.')[1]) for g in hfbad]
                    badnets = [str(g.split('.')[3].split(' - ')[0]) for g in hfbad]
                    goodstas = [str(trace.stats.station) for trace in st_hfgood]
                    goodchans = [str(trace.stats.channel) for trace in st_hfgood]
                    goodnets = [str(trace.stats.network) for trace in st_hfgood]
            elif feedback.lower() == 'y':
                # put info in database
                badstas = [str(g.split('.')[0]) for g in hfbad]
                badchans = [str(g.split('.')[1]) for g in hfbad]
                badnets = [str(g.split('.')[3].split(' - ')[0]) for g in hfbad]
                goodstas = [str(trace.stats.station) for trace in st_hfgood]
                goodchans = [str(trace.stats.channel) for trace in st_hfgood]
                goodnets = [str(trace.stats.network) for trace in st_hfgood]
                maxreachedHF = True
                # Insert maxdisthf into database
                with connection:
                    cursor = connection.cursor()
                    try:
                        cursor.execute('UPDATE events SET maxdistHF_km = ? WHERE Eid = ?', (maxdist, event_id))
                    except Exception as e:
                        print(e)
            # Fill info into database
            Sids_bad = []
            Sids_good = []
            with connection:
                cursor = connection.cursor()
                for i, badsta in enumerate(badstas):
                    try:
                        cursor_output = cursor.execute('SELECT Sid FROM stations WHERE Name=? AND Network=? AND Channel=?',
                                                       (badstas[i], badnets[i], badchans[i]))
                        retrieved_data = cursor_output.fetchall()
                        if retrieved_data:
                            Sids_bad += [val[0] for val in retrieved_data]
                    except Exception as e:
                        print(e)
                for i, goodsta in enumerate(goodstas):
                    try:
                        cursor_output = cursor.execute('SELECT Sid FROM stations WHERE Name=? AND Network=? AND Channel=?',
                                                       (goodstas[i], goodnets[i], goodchans[i]))
                        retrieved_data = cursor_output.fetchall()
                        if retrieved_data:
                            Sids_good += [val[0] for val in retrieved_data]
                    except Exception as e:
                        print(e)
                if Sids_bad:
                    for Sid in Sids_bad:
                        cursor.execute('UPDATE sta_nearby SET detect_HF = ? WHERE station_id = ? AND event_id = ?',
                                       (0, Sid, event_id))
                if Sids_good:
                    for Sid in Sids_good:
                        cursor.execute('UPDATE sta_nearby SET detect_HF = ? WHERE station_id = ? AND event_id = ?',
                                       (1, Sid, event_id))
                if maxreachedHF is True:
                    cursor.execute('UPDATE sta_nearby SET detect_HF = ? WHERE stasource_radius_km > ? AND event_id = ?',
                                   (0, float(maxdist), event_id))

        if maxreachedLP is False:
            st_lp = st.copy()
            st_lp = st_lp.select(channel='B*') + st_lp.select(channel='H*') + st_lp.select(channel='L*')
            #reorder by distance
            st_lp = st_lp.sort(keys=['rdist', 'channel'])
            temp = Stream()
            cosfilt = (0.5/LPlims[1], 1/LPlims[1], 1/LPlims[0], 2/LPlims[0])
            for i, trace in enumerate(st_lp):
                try:
                    trace.remove_response(output=LPoutput, pre_filt=cosfilt)
                    temp = temp + trace
                except:
                    print('Failed to remove response for %s, deleting this station' % (trace.stats.station + trace.stats.channel,))
            st_lp = temp

            print ('Interactive review of long periods starting\n')
            # Interactive review of long periods
            zp = reviewData.InteractivePlot(st_lp, maxtraces=maxtraces, textline=['Station corrected LP data for your review',
                                            'Delete bad traces and take note of distance when you reach maximum observation',
                                            'hit Q to exit when ready to continue'])
            # Get info about what was kept and what wasn't
            st_lpgood = zp.st_current
            lpbad = zp.deleted

            connection = None
            connection = lite.connect(database)
            maxdist = np.max([trace.stats.rdist for trace in st_lpgood])
            feedback = input(('is %4.1f km the maximum distance observed for lp? Y or N\n') % (maxdist,))
            if feedback.lower() == 'n':
                feedback = input('was the maximum distance less than %4.1f km?\n' % (maxdist,))
                if feedback.lower() == 'y':
                    maxdist = float(input('Enter the maximum distance observed in km (round up to nearest km)\n'))
                    # delete all stations that were greater than maximum distance and that are in st_lpbad, and put good stations in database
                    goodstas = [str(trace.stats.station) for trace in st_lpgood if trace.stats.rdist < maxdist]
                    goodchans = [str(trace.stats.channel) for trace in st_lpgood if trace.stats.rdist < maxdist]
                    goodnets = [str(trace.stats.network) for trace in st_lpgood if trace.stats.rdist < maxdist]
                    badstas = [str(trace.stats.station) for trace in st_lpgood if trace.stats.rdist > maxdist] + [str(g.split('.')[0]) for g in lpbad]
                    badchans = [str(trace.stats.channel) for trace in st_lpgood if trace.stats.rdist > maxdist] + [str(g.split('.')[1]) for g in lpbad]
                    badnets = [str(trace.stats.network) for trace in st_lpgood if trace.stats.rdist > maxdist] + [str(g.split('.')[3].split(' - ')[0]) for g in lpbad]
                    maxreachedLP = True
                    # Insert maxdistlp into database
                    with connection:
                        cursor = connection.cursor()
                        try:
                            cursor.execute('UPDATE events SET maxdistLP_km = ? WHERE Eid = ?', (maxdist, event_id))
                        except Exception as e:
                            print(e)
                if feedback.lower() == 'n':
                    print('Filling current info into database, then will load more distant data for further analysis\n')
                    badstas = [str(g.split('.')[0]) for g in lpbad]
                    badchans = [str(g.split('.')[1]) for g in lpbad]
                    badnets = [str(g.split('.')[3].split(' - ')[0]) for g in lpbad]
                    goodstas = [str(trace.stats.station) for trace in st_lpgood]
                    goodchans = [str(trace.stats.channel) for trace in st_lpgood]
                    goodnets = [str(trace.stats.network) for trace in st_lpgood]
            elif feedback.lower() == 'y':
                # put info in database
                badstas = [str(g.split('.')[0]) for g in lpbad]
                badchans = [str(g.split('.')[1]) for g in lpbad]
                badnets = [str(g.split('.')[3].split(' - ')[0]) for g in lpbad]
                goodstas = [str(trace.stats.station) for trace in st_lpgood]
                goodchans = [str(trace.stats.channel) for trace in st_lpgood]
                goodnets = [str(trace.stats.network) for trace in st_lpgood]
                maxreachedLP = True
                # Insert maxdistlp into database
                with connection:
                    cursor = connection.cursor()
                    try:
                        cursor.execute('UPDATE events SET maxdistLP_km = ? WHERE Eid = ?', (maxdist, event_id))
                    except Exception as e:
                        print(e)
            # Fill info into database
            Sids_bad = []
            Sids_good = []
            with connection:
                cursor = connection.cursor()
                for i, badsta in enumerate(badstas):
                    try:
                        cursor_output = cursor.execute('SELECT Sid FROM stations WHERE Name=? AND Network=? AND Channel=?',
                                                       (badstas[i], badnets[i], badchans[i]))
                        retrieved_data = cursor_output.fetchall()
                        if retrieved_data:
                            Sids_bad += [val[0] for val in retrieved_data]
                    except Exception as e:
                        print(e)
                for i, goodsta in enumerate(goodstas):
                    try:
                        cursor_output = cursor.execute('SELECT Sid FROM stations WHERE Name=? AND Network=? AND Channel=?',
                                                       (goodstas[i], goodnets[i], goodchans[i]))
                        retrieved_data = cursor_output.fetchall()
                        if retrieved_data:
                            Sids_good += [val[0] for val in retrieved_data]
                    except Exception as e:
                        print(e)
                if Sids_bad:
                    for Sid in Sids_bad:
                        cursor.execute('UPDATE sta_nearby SET detect_LP = ? WHERE station_id = ? AND event_id = ?',
                                       (0, Sid, event_id))
                if Sids_good:
                    for Sid in Sids_good:
                        cursor.execute('UPDATE sta_nearby SET detect_LP = ? WHERE station_id = ? AND event_id = ?',
                                       (1, Sid, event_id))
                if maxreachedLP is True:
                    cursor.execute('UPDATE sta_nearby SET detect_LP = ? WHERE stasource_radius_km > ? AND event_id = ?',
                                   (0, float(maxdist), event_id))


def make_measurementsHF(event_id, buffer_sec=100., HFlims=(1., 5.), HFoutput='VEL',
                        minradius=0., maxradius=None, maxtraces=15, taper=None, detrend='demean',
                        database='/Users/kallstadt/LSseis/landslideDatabase/lsseis.db',
                        path='/Users/kallstadt/LSseis/landslideDatabase'):
    """
    Make measurements on confirmed high frequency (HF) detections, pulls data from IRIS and other sources
    and links station correction info, then displays in Interactive Plot. User should then delete all traces
    that are clipped, then press C to perform station correction, Make amplitude picks at start and
    end of each signal (best guess), then quit interactive plotting
    by pressing q and and it will then use the time of the amplitude picks to
    populate database with several measurements about each signal. Do not make
    measurements on wonky signals.
    Note: This must be done in ipython with pylab enabled for interactive plotting to work. May also require
    that an interactive backend be used for matplotlib (Qt5Agg works).

    Args:
        event_id (int): id of event to review
        buffer_sec (float): Number of seconds to add on either end of start and end time (makes viewing easier)
        HFlims: tuple or list of lower and upper frequency limits in Hz for HF filtering (1-5 Hz standard)
        HFoutput (str): Output type for HF station correction ('VEL' standard, could also be 'DISP' or 'ACC')
        minradius (float): minimum distance in km to search for stations (from source)
        maxradius (float): maximum distance in km to search for stations (from source)
        maxtraces (int): Number of traces to view at a time
        taper (0.05): percentage taper (or None for no tapering)
        detrend (str): type of detrending to use. 'demean' or 'linear' are common options
        database (str): Full file path of database file location
        path (str): Path to location of sac files (upstream from relative file paths listed in database)

    """
    # built cosine filter that will be used
    cosfilt = (0.5*HFlims[0], HFlims[0], HFlims[1], 2*HFlims[1])

    # Get event info
    evDict = findsta.getEventInfo(event_id, database=database)
    print((('Now analysing Eid %s - %s') % (event_id, evDict['Name'])))

    # Find stations where detect_HF=1 for this event
    staDict = findsta.getStaInfo(event_id, database=database, detectHF=True, minradius=minradius, maxradius=maxradius)
    #datlocs = reviewData.unique_list([staDict[k]['source'] for k in staDict])
    sttemp = Stream()
    if 'IRIS' in evDict['DatLocation']:
        stalist = []
        for sdct in staDict:
            if 'IRIS' in sdct['source']:
                if 'AlaskaData' not in evDict['DatLocation']:
                    stalist.append((sdct['Name'], sdct['Channel'], sdct['Network'], '*'))
                else:
                    if 'AV' not in sdct['Network'] and 'AK' not in sdct['Network']:
                        stalist.append((sdct['Name'], sdct['Channel'], sdct['Network'], '*'))
        if len(stalist) != 0:
            try:
                sttemp += reviewData.getdata_exact(stalist, evDict['StartTime'] - buffer_sec, evDict['EndTime'] + buffer_sec,
                                                   attach_response=True, clientname='IRIS', detrend='demean')
            except Exception as e:
                print(e)
    if 'NCEDC' in evDict['DatLocation']:
        stalist = [(sdct['Name'], sdct['Channel'], sdct['Network'], '*') for sdct in staDict if 'NCEDC' in sdct['source']]
        if len(stalist) != 0:
            sttemp += reviewData.getdata_exact(stalist, evDict['StartTime'] - buffer_sec, evDict['EndTime'] + buffer_sec,
                                               attach_response=True, clientname='NCEDC', detrend='demean')
    if evDict['DatLocation'] is None:
        print('You need to populate the DatLocation field for this event, no sac files loaded')
    if 'sac' in evDict['DatLocation']:
        datloc1 = evDict['DatLocation'].split(',')
        datloc1 = [x.strip() for x in datloc1 if 'sac' in x]
        for datl in datloc1:
            fullpath = os.path.join(path, datl.split(':')[1])
            filenames = glob.glob(fullpath)
            if len(filenames) > 0:
                if 'AlaskaData' in datl:  # Need to do some cheating to attach response info if Iliamna sac data - attach oldest response info available at IRIS for each station
                    # Only keep filenames of stations we want to read in
                    newfilenames = []
                    namelist = [sdct['Name'] for sdct in staDict]
                    for filen in filenames:
                        nam = filen.split('/')[-1].split('.')[0]
                        if nam in namelist:
                            newfilenames.append(filen)
                    stsac = reviewData.getdata_sac(newfilenames, attach_response=True, starttime=evDict['StartTime']-buffer_sec,
                                                   endtime=evDict['EndTime']+buffer_sec, detrend='demean')
                    stsac = Stream([trace for trace in stsac if trace.max() != 0.0])  # Get rid of any empty ones
                    originalstt = []
                    # Get earliest start time from IRIS for AK and AV
                    url1 = ('http://service.iris.edu/fdsnws/station/1/query?network=AV&level=station&format=text&nodata=404')
                    url2 = ('http://service.iris.edu/fdsnws/station/1/query?network=AK&level=station&format=text&nodata=404')
                    f = urllib.request.urlopen(url1)
                    file1 = f.read()
                    lines1 = [line.split('|') for line in file1.split('\n')[1:]]
                    stdict = {}
                    f.close()
                    f = urllib.request.urlopen(url2)
                    file2 = f.read()
                    lines2 = [line.split('|') for line in file2.split('\n')[1:]]
                    f.close()
                    for line in lines1[:-1]:
                        stdict[line[1]] = (line[0], UTCDateTime(line[6]))
                    for line in lines2[:-1]:
                        if line[1] in list(stdict.keys()):
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
                            print('could not attach response info for %s, station correction will not work' % temp.stats.station)
                            continue
                        if 'response' not in list(temp.stats.keys()):
                            # Change starttime of temp to earliest start time at IRIS for this station
                            try:
                                temp.stats.starttime = stdict[temp.stats.station][1] + 86400.
                            except:
                                print('could not attach response info for %s, station correction will not work' % temp.stats.station)
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
                                                   endtime=evDict['EndTime']+buffer_sec, detrend='demean')
                sttemp += stsac

                for trace in sttemp:
                    if '--' in trace.stats.location:
                        trace.stats.location = ''

    # Delete any that are not in list
    nets, stas, chans = list(zip(*[[sdct['Network'], sdct['Name'], sdct['Channel']] for sdct in staDict]))  # Exclude location code because inconsistencies here can cause data to not be found
    st = Stream()
    for n1, s1, c1 in zip(nets, stas, chans):
        st += sttemp.select(station=s1, network=n1, channel=c1)
    if len(st) == 0 and event_id == 106:
        for s1, c1 in zip(stas, chans):
            st += sttemp.select(station=s1, channel=c1)
        st.merge()
    # Attach distaz
    st = findsta.attach_distaz(st, evDict['Latitude'], evDict['Longitude'], database=database)
    st = st.sort(keys=['rdist', 'channel'])

    # Preprocess
    if detrend is not None:
        st.detrend(detrend)
    if taper is not None:
        st.taper(max_percentage=taper, type='cosine')

    # Interactive plotting to pick start and end time and make amplitude picks
    textline = ['Delete any that are clipped after making amplitude picks on them to get start and end times',
                'Then press c for station corrected VEL data', 'Make amplitude picks by hitting A at start and end of signal',
                'these times will be used for other calculations', 'Zoom out all the way before exiting']
    zp = reviewData.InteractivePlot(st, cosfilt=cosfilt, output=HFoutput, textline=textline, maxtraces=maxtraces, taper=taper)

    # Get amplitude pick times
    picks = zp.picks
    clipped = zp.deleted
    stcorr = zp.st_current

    # Make sure only using amplitude picks
    amppicks = [picks[key] for key in picks if picks[key]['phase'] in 'A' if picks[key]['stachan'] not in clipped]

    # Separate amplitude picks that are valid from those that were clipped (deleted)
    durpickonly = [picks[key] for key in picks if picks[key]['phase'] in 'A' if picks[key]['stachan'] in clipped]

    # Use start and end time from amplitude pick to compute the other things
    starttimes = []
    endtimes = []
    p2pamp = []
    astas = []
    anets = []
    achans = []
    alocs = []
    aamps = []
    durations = []
    sliceorig = Stream()
    presliceorig = Stream()
    for amppick in amppicks:
        starttimes.append(amppick['picktime'][0])
        endtimes.append(amppick['picktime'][1])
        durations.append(endtimes[-1] - starttimes[-1])
        code = amppick['stachan'].split('-')[0].strip().split('.')
        code1 = '.'.join((code[3], code[0], code[2], code[1]))
        astas.append(code[0])
        anets.append(code[3])
        achans.append(code[1])
        alocs.append(code[2])
        aamps.append(amppick['weight'])
        temp1 = stcorr.select(id=str(code1)).copy()
        #stt = temp1[0].stats.starttime
        temp = temp1.slice(amppick['picktime'][0], amppick['picktime'][1])
        p2pamp.append(np.max(temp[0].data) - np.min(temp[0].data))
        sliceorig += st.select(id=str(code1)).copy().slice(amppick['picktime'][0], amppick['picktime'][1])
        presliceorig += st.select(id=str(code1)).copy().slice(st.select(id=str(code1))[0].stats.starttime, amppick['picktime'][0])
    peakfreq = sigproc.peakfreq(sliceorig)
    meansqfreq, var = sigproc.meansqfreq(sliceorig)
    meansqfreqSN, var = sigproc.meansqfreqSN(sliceorig, presliceorig, SNrat=2.0)
    #domfreq = sigproc.domfreq(sliceorig)

    stations, networks, channels, locations, SRid = list(zip(*[[sdct['Name'], sdct['Network'], sdct['Channel'],
                                                        sdct['LocationCode'], sdct['SRid']] for sdct in staDict]))

    # Put in the database
    connection = None
    connection = lite.connect(database)
    with connection:
        cursor = connection.cursor()
        for i in np.arange(len(astas)):
            # get SRid
            thisSRid = [sid for j, sid in enumerate(SRid) if astas[i] in stations[j] and achans[i] in channels[j] and alocs[i] in locations[j] and anets[i] in networks[j]]
            if len(thisSRid) == 0:
                thisSRid = [sid for j, sid in enumerate(SRid) if astas[i] in stations[j] and achans[i] in channels[j] and anets[i] in networks[j]]
            try:
                if len(thisSRid) > 1:
                    station5 = '%s.%s.%s.%s' % (astas[i], achans[i], alocs[i], anets[i])
                    print(('Found multiple matches when excluding location code for %s, entering data into all of them' % station5))
                for thisid in thisSRid:
                    cursor.execute('UPDATE sta_nearby SET starttimeHF = ?, endtimeHF = ?, peakfreqraw = ?, meansqfreqraw = ?, meansqfreqSNRraw = ?, duration_secHF = ?, absmaxampHF = ?, p2pmaxampHF = ? WHERE SRid = ?', (starttimes[i].strftime('%Y-%m-%d %H:%M:%S'), endtimes[i].strftime('%Y-%m-%d %H:%M:%S'), peakfreq[i], meansqfreq[i], meansqfreqSN[i], durations[i], aamps[i], p2pamp[i], thisid))
            except Exception as e:
                print(e)
                print(('Could not put %s in database' % station5))

    # Put duration measurements in the database too
    dstarttimes = []
    dendtimes = []
    ddurations = []
    dstas = []
    dnets = []
    dchans = []
    dlocs = []
    for durpick in durpickonly:
        dstarttimes.append(durpick['picktime'][0])
        dendtimes.append(durpick['picktime'][1])
        ddurations.append(dendtimes[-1] - dstarttimes[-1])
        code = durpick['stachan'].split('-')[0].strip().split('.')
        code1 = '.'.join((code[3], code[0], code[2], code[1]))
        dstas.append(code[0])
        dnets.append(code[3])
        dchans.append(code[1])
        dlocs.append(code[2])
    connection = None
    connection = lite.connect(database)
    with connection:
        cursor = connection.cursor()
        for i in np.arange(len(dstas)):
            # get SRid
            thisSRid = [sid for j, sid in enumerate(SRid) if dstas[i] in stations[j] and dchans[i] in channels[j] and dlocs[i] in locations[j] and dnets[i] in networks[j]]
            try:
                cursor.execute('UPDATE sta_nearby SET starttimeHF = ?, endtimeHF = ?, duration_secHF = ? WHERE SRid = ?', (dstarttimes[i].strftime('%Y-%m-%d %H:%M:%S'), dendtimes[i].strftime('%Y-%m-%d %H:%M:%S'), ddurations[i], thisSRid[0]))
            except Exception as e:
                print(e)


def make_measurementsLP(event_id, buffer_sec=100., LPlims=(20., 60.), LPoutput='DISP', taper=0.05, detrend='linear',
                        minradius=0., maxradius=None, maxtraces=15, path='/Users/kallstadt/LSseis/landslideDatabase',
                        database='/Users/kallstadt/LSseis/landslideDatabase/lsseis.db'):
    """
    Make measurements on confirmed long period detections, pulls data from IRIS and other sources
    and links station correction info, then displays in Interactive Plot. User should then delete all traces
    that are clipped, then press C to perform station correction, Make amplitude picks at start and
    end of each signal (best guess), then quit interactive plotting
    by pressing q and and it will then use the time of the amplitude picks to
    populate database with several measurements about each signal. Do not make
    measurements on wonky signals.
    Note: This must be done in ipython with pylab enabled for interactive plotting to work. May also require
    that an interactive backend be used for matplotlib (Qt5Agg works).

    Args:
        event_id (int): id of event to review
        buffer_sec (float): Number of seconds to add on either end of start and end time (makes viewing easier)
        LPlims: tuple or list of lower and upper frequency limits in seconds for HF filtering (20-60 sec standard)
        LPoutput (str): Output type for HF station correction ('DISP' standard, could also be 'VEL' or 'ACC')
        taper (0.05): percentage taper (or None for no tapering)
        detrend (str): type of detrending to use. 'demean' or 'linear' are common options
        minradius (float): minimum distance in km to search for stations (from source)
        maxradius (float): maximum distance in km to search for stations (from source)
        maxtraces (int): Number of traces to view at a time
        path (str): Path to location of sac files (upstream from relative file paths listed in database)
        database (str): Full file path of database file location

    """
    # built cosine filter that will be used
    cosfilt = (0.5/LPlims[1], 1/LPlims[1], 1/LPlims[0], 2/LPlims[0])

    # Get event info
    evDict = findsta.getEventInfo(event_id, database=database)
    print((('Now analysing Eid %s - %s') % (event_id, evDict['Name'])))

    # Find stations where detect_HF=1 for this event
    staDict = findsta.getStaInfo(event_id, database=database, detectLP=True, minradius=minradius, maxradius=maxradius)
    datlocs = reviewData.unique_list([sdct['source'] for sdct in staDict])
    sttemp = Stream()
    if 'IRIS' in evDict['DatLocation'] or 'IRIS' in datlocs:
        stalist = []
        for sdct in staDict:
            if 'IRIS' in sdct['source']:
                if 'AlaskaData' not in evDict['DatLocation']:
                    stalist.append((sdct['Name'], sdct['Channel'], sdct['Network'], '*'))
                else:
                    if 'AV' not in sdct['Network'] and 'AK' not in sdct['Network']:
                        stalist.append((sdct['Name'], sdct['Channel'], sdct['Network'], '*'))
        if len(stalist) != 0:
            try:
                sttemp += reviewData.getdata_exact(stalist, evDict['StartTime'] - buffer_sec, evDict['EndTime'] + buffer_sec,
                                                   attach_response=True, clientname='IRIS', detrend=detrend)
            except Exception as e:
                print(e)
    if 'NCEDC' in evDict['DatLocation'] or 'NCEDC' in datlocs:
        stalist = [(sdct['Name'], sdct['Channel'], sdct['Network'], '*') for sdct in staDict if 'NCEDC' in sdct['source']]
        if len(stalist) != 0:
            sttemp += reviewData.getdata_exact(stalist, evDict['StartTime'] - buffer_sec, evDict['EndTime'] + buffer_sec,
                                               attach_response=True, clientname='NCEDC', detrend=detrend)
    if evDict['DatLocation'] is None:
        print('You need to populate the DatLocation field for this event, no sac files loaded')
    if 'sac' in evDict['DatLocation']:
        datloc1 = evDict['DatLocation'].split(',')
        datloc1 = [x.strip() for x in datloc1 if 'sac' in x]
        for datl in datloc1:
            fullpath = os.path.join(path, datl.split(':')[1])
            filenames = glob.glob(fullpath)
            if len(filenames) > 0:
                if 'AlaskaData' in datl:  # Need to do some cheating to attach response info if Iliamna sac data - attach oldest response info available at IRIS for each station
                    # Only keep filenames of stations we want to read in
                    newfilenames = []
                    namelist = [sdct['Name'] for sdct in staDict]
                    for filen in filenames:
                        nam = filen.split('/')[-1].split('.')[0]
                        if nam in namelist:
                            newfilenames.append(filen)
                    stsac = reviewData.getdata_sac(newfilenames, attach_response=True, starttime=evDict['StartTime']-buffer_sec,
                                                   endtime=evDict['EndTime']+buffer_sec, detrend=detrend)
                    stsac = Stream([trace for trace in stsac if trace.max() != 0.0])  # Get rid of any empty ones
                    originalstt = []
                    # Get earliest start time from IRIS for AK and AV
                    url1 = ('http://service.iris.edu/fdsnws/station/1/query?network=AV&level=station&format=text&nodata=404')
                    url2 = ('http://service.iris.edu/fdsnws/station/1/query?network=AK&level=station&format=text&nodata=404')
                    f = urllib.request.urlopen(url1)
                    file1 = f.read()
                    lines1 = [line.split('|') for line in file1.split('\n')[1:]]
                    stdict = {}
                    f.close()
                    f = urllib.request.urlopen(url2)
                    file2 = f.read()
                    lines2 = [line.split('|') for line in file2.split('\n')[1:]]
                    f.close()
                    for line in lines1[:-1]:
                        stdict[line[1]] = (line[0], UTCDateTime(line[6]))
                    for line in lines2[:-1]:
                        if line[1] in list(stdict.keys()):
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
                            print('could not attach response info for %s, station correction will not work' % temp.stats.station)
                            continue
                        if 'response' not in list(temp.stats.keys()):
                            # Change starttime of temp to earliest start time at IRIS for this station
                            try:
                                temp.stats.starttime = stdict[temp.stats.station][1] + 86400.
                            except:
                                print('could not attach response info for %s, station correction will not work' % temp.stats.station)
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
                                                   endtime=evDict['EndTime']+buffer_sec, detrend='linear')
                sttemp += stsac

    for trace in sttemp:
        if '--' in trace.stats.location:
            trace.stats.location = ''

    # Delete any that are not in list
    nets, stas, chans = list(zip(*[[sdct['Network'], sdct['Name'], sdct['Channel']] for sdct in staDict]))  # Exclude location code because inconsistencies here can cause data to not be found
    st = Stream()
    for n1, s1, c1 in zip(nets, stas, chans):
        st += sttemp.select(station=s1, network=n1, channel=c1)
    # Attach distaz
    st = findsta.attach_distaz(st, evDict['Latitude'], evDict['Longitude'], database=database)
    st = st.sort(keys=['rdist', 'channel'])

    # Attach distaz
    st = findsta.attach_distaz(st, evDict['Latitude'], evDict['Longitude'], database=database)
    st = st.sort(keys=['rdist', 'channel'])

    # Preprocess (already detrended, if option was chosen, when downloaded data)
    if taper is not None:
        st.taper(max_percentage=taper, type='cosine')

    # Perform station correction
    temp = Stream()
    for i, trace in enumerate(st):
        try:
            trace.remove_response(output=LPoutput, pre_filt=cosfilt)
            temp = temp + trace
        except:
            print('Failed to remove response for %s, deleting this station' % (trace.stats.station + trace.stats.channel,))
    st = temp

    # Interactive plotting to pick start and end time and make amplitude picks
    textline = ['Make amplitude picks on good signals, check and make sure amplitudes are not crazy']
    zp = reviewData.InteractivePlot(st, textline=textline, taper=taper)

    # Get amplitude pick times
    picks = zp.picks
    clipped = zp.deleted

    # Make sure only using amplitude picks
    amppicks = [picks[key] for key in picks if picks[key]['phase'] in 'A' if picks[key]['stachan'] not in clipped]

    # Use start and end time from amplitude pick to compute the other things
    starttimes = []
    endtimes = []
    p2pamp = []
    astas = []
    anets = []
    achans = []
    alocs = []
    aamps = []
    durations = []
    sliceorig = Stream()
    presliceorig = Stream()
    for amppick in amppicks:
        starttimes.append(amppick['picktime'][0])
        endtimes.append(amppick['picktime'][1])
        durations.append(endtimes[-1] - starttimes[-1])
        code = amppick['stachan'].split('-')[0].strip().split('.')
        code1 = '.'.join((code[3], code[0], code[2], code[1]))
        astas.append(code[0])
        anets.append(code[3])
        achans.append(code[1])
        alocs.append(code[2])
        aamps.append(amppick['weight'])
        temp1 = st.select(id=str(code1)).copy()
        temp = temp1.slice(amppick['picktime'][0], amppick['picktime'][1])
        p2pamp.append(np.max(temp[0].data) - np.min(temp[0].data))
        sliceorig += st.select(id=str(code1)).copy().slice(amppick['picktime'][0], amppick['picktime'][1])
        presliceorig += st.select(id=str(code1)).copy().slice(st.select(id=str(code1))[0].stats.starttime, amppick['picktime'][0])

    stations, networks, channels, locations, SRid = list(zip(*[[sdct['Name'], sdct['Network'], sdct['Channel'],
                                                         sdct['LocationCode'], sdct['SRid']] for sdct in staDict]))

    # Put in the database
    connection = None
    connection = lite.connect(database)
    with connection:
        cursor = connection.cursor()
        for i in np.arange(len(astas)):
            # get SRid
            thisSRid = [sid for j, sid in enumerate(SRid) if astas[i] in stations[j] and achans[i] in channels[j] and alocs[i] in locations[j] and anets[i] in networks[j]]
            if len(thisSRid) == 0:
                thisSRid = [sid for j, sid in enumerate(SRid) if astas[i] in stations[j] and achans[i] in channels[j] and anets[i] in networks[j]]
            try:
                if len(thisSRid) > 1:
                    station5 = '%s.%s.%s.%s' % (astas[i], achans[i], alocs[i], anets[i])
                    print(('Found multiple matches when excluding location code for %s, entering data into all of them' % station5))
                for thisid in thisSRid:
                    cursor.execute('UPDATE sta_nearby SET starttimeLP = ?, endtimeLP = ?, duration_secLP = ?, absmaxampLP = ?, p2pmaxampLP = ? WHERE SRid = ?', (starttimes[i].strftime('%Y-%m-%d %H:%M:%S'), endtimes[i].strftime('%Y-%m-%d %H:%M:%S'), durations[i], aamps[i], p2pamp[i], thisid))
            except Exception as e:
                print(e)
                print(('Could not put %s in database' % station5))


def populate_redoubt_stanearby(event_ids, database='/Users/kallstadt/LSseis/landslideDatabase/lsseis.db',
                               update=True):
    """This function is only for a specific set of lahar events that occurred during the Redoubt eruption
    """
    if type(event_ids) is int:
        event_ids = [event_ids]
    for event_id in event_ids:
        evdict = findsta.getEventInfo(event_id, database=database)

    stations = ['RD01', 'RD02', 'RD03', 'RDW']

    connection = None
    connection = lite.connect(database)
    cursor = connection.cursor()

    #get Sid and lat lon for each entry in line
    val = 0
    val1 = 0
    val2 = 0
    valbad = 0
    for event_id in event_ids:
        for sta in stations:
                #get the Sid from stations table
                cursor_output = cursor.execute('SELECT Sid,Latitude,Longitude FROM stations WHERE Name=?', (sta,))
                retrieved_data = cursor_output.fetchall()
                for dat in retrieved_data:
                    Sid, sta_lat, sta_lon = dat
                    try:
                        cursor_output = cursor.execute('SELECT SRid FROM sta_nearby WHERE event_id=? AND station_id=?',
                                                       (event_id, Sid))
                        temp = cursor_output.fetchall()
                    except:
                        temp = []
                    if len(temp) == 0:
                        #if it isn't already there, use iris's distaz webservice to calculate stuff and save it in the table
                        #result = client.distaz(sta_lat, sta_lon, event_lat, event_lon)
                        backazimuth, azimuth, distance = reviewData.pyproj_distaz(sta_lat, sta_lon, evdict['Latitude'], evdict['Longitude'])
                        with connection:
                            try:
                                cursor.execute('INSERT INTO sta_nearby(event_id,station_id,stasource_radius_km,az,baz) VALUES(?,?,?,?,?)',
                                               (event_id, Sid, '%3.3f' % distance, '%3.2f' % azimuth, '%3.2f' % backazimuth))
                                val += 1
                            except Exception as f:
                                print(f)
                                valbad += 1
                    elif update is True:
                        # update entry, if it is there
                        with connection:
                            try:
                                #result = client.distaz(sta_lat, sta_lon, event_lat, event_lon)
                                backazimuth, azimuth, distance = reviewData.pyproj_distaz(sta_lat, sta_lon, evdict['Latitude'], evdict['Longitude'])
                                cursor.execute("UPDATE sta_nearby SET stasource_radius_km = ?, az = ?, baz = ? WHERE event_id =? AND station_id =?;",
                                               ('%3.3f' % distance, '%3.2f' % azimuth, '%3.2f' % backazimuth, event_id, Sid))
                                val1 += 1
                            except Exception as f:
                                print(f)
                                valbad += 1
                    else:
                        val2 += 1

        print(('added %s entries to sta_nearby table, updated %s entries, %s left as is, %s not added because of error' % (val, val1, val2, valbad)))
