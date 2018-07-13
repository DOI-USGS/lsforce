#!/usr/bin/env python

"""
Functions for interacting with and populating the seismogenic landslide database
"""

import sqlite3 as lite
from obspy import Stream, UTCDateTime
import numpy as np
from reviewData import reviewData


def getEventInfo(event_id, database='/Users/kallstadt/LSseis/landslideDatabase/lsseis.db'):
    """This function pulls the entire entry from the events table into a dictionary (e.g. Lats/Lons, volumes, everything in row)

    Args:

        event_id (int): id number corresponding to the event of interest in the database
        database (str): Full file path to SQLite3 landslide database file

    Returns:
       Dictionary containing all entries of all fields from event table of database

    Dictionary Fields:
    'Eid' - event id
    'Name' - working name of event (non-unique)
    'StartTime' - starttime as UTCDateTime
    'EndTime' - endtime as UTCDateTime
    'Latitude' - Latitude of approximate center of landslide
    'Longitude'
    'LocUncert_km' - approximate location uncertainty in km
    'Top_lat' - Latitude of top (highest point)
    'Top_lon'
    'Toe_lat' - Latitude of toe (furthest point)
    'Toe_lon'
    'Type' - Type of event
    'Area_source' - Best estimate of area of the source area in m**2
    'Area_source_low' - Lower bound estimate of area of the source area in m**2
    'Area_source_high'- Higher bound of area of the source area in m**2
    'Volume' - Best estimate of volume in m**3
    'Volume_low' - Lower bound estimate of volume in m**3
    'Volume_high' - Higher bound estimate of volume in m**3
    'Mass' - Best estimate of mass in kg
    'Mass_low' - Lower bound estimate of mass in kg
    'Mass_high' - Higher bound estimate of mass in kg
    'Area_total' - Best estimate of total area in m**2
    'Area_total_low' - Lower bound estimate of total area in m**2
    'Area_total_high' - Higher bound estimate of total area in m**2
    'H' - Best estimate drop height from top to toe in meters
    'H_low' - Lower bound estimate drop height from top to toe in meters
    'H_high' - Upper bound estimate drop height from top to toe in meters
    'L' - Best estimate horizontal travel distance from top to toe in meters
    'L_low' - Lower bound estimate horizontal travel distance from top to toe in meters
    'L_high' - Upper bound estimate horizontal travel distance from top to toe in meters
    'OtherDataQuality1to5' - Subjective rating of ancillary data (volume etc.), 1 low to 5 high
    'LPpotential' - Potential for long period analysis, 0 = No potential, 1 = long period waves observable and usable, 2 = long periods waves observable but not usable
    'maxdistHF_km' - Maximum distance high frequency (1-5 Hz, VEL) seismic waves are observed in km
    'maxdistLP_km' - Maximum distance long period (20-60 sec, DISP) seismic waves are observed in km
    'A0' - MAY DELETE THE NEXT FIVE, DIDNT WORK WELL, FROM FIT TO ATTENUATION
    'A01std' -
    'atten' -
    'atten1std' -
    'fmedian' - Median frequency (Hz)
    'DatLocation' - Location of seismic data, IRIS, sac files with relative file path
    'Seismic_observed' - 1 if seismic signals found, 0 if not.
    """
    connection = None
    connection = lite.connect(database, detect_types=lite.PARSE_DECLTYPES)
    connection.text_factory = str
    with connection:
        cursor = connection.cursor()
        try:
            event_id = int(event_id)
            cursor_output = (cursor.execute('SELECT * FROM events WHERE Eid = ?', (event_id, )))
            dat = cursor_output.fetchall()[0]
            temp = cursor.description
            names = [t[0] for t in temp]
            eventDict = dict(list(zip(names, dat)))
            # Convert to the right units
            temp = eventDict['StartTime'].split(' ')
            eventDict['StartTime'] = UTCDateTime(temp[0]+'T'+temp[1])
            temp = eventDict['EndTime'].split(' ')
            eventDict['EndTime'] = UTCDateTime(temp[0]+'T'+temp[1])
        except Exception as e:
            print(e)
            return

        return eventDict


def getStaInfo(event_id, maxradius=None, minradius=0., detectHF=None, detectLP=None, both=False,
               ampdetectHF=None, ampdetectLP=None, durationHF=None, durationLP=None,
               chanuse='*', numstas=None, database='/Users/kallstadt/LSseis/landslideDatabase/lsseis.db'):
    """Get all the station detection information with or without constraints on radius,
    whether or not it was detected in HF or LP bands, and whether or not there
    are amplitude and duration measurements in HF or LP bands.

    Args:
        event_id (int): id number corresponding to the event of interest in the database
        maxradius (float): Maximum radius to search in km
        detectHF: Return HF detections? None, True, or False
        detectLP: Return LP detections? None, True, or False
        both (bool): True if want both detectHF and detectLP, but don't need both to be true for a given station
        ampdetectHF: Return stations with HF amplitude measurements? None, True, or False
        ampdetectLP: Return stations with HF amplitude measurements? None, True, or False
        durationHF: Return stations with HF duration measurements? None, True or False
        durationLP: Return stations with LP duration measurements? None, True or False
        chanuse: Single string of channel names, return only stations that match channels defined,
            '*' for all channels
        numstas (float): Number of detections meeting criteria to return, starting from closest distance,
            if None, will return all. Channels from the same station count as one.
        database (str): Full file path to SQLite3 landslide database file

    Returns:
        List of dictionaries containing all entries of all fields from event table
            of database

        Dictionary entries
        'Sid' - Station id
        'Name' - Station name
        'Network' - Network code
        'Channel' - Channel code
        'LocationCode' - Location code
        'Latitude' - Station latitude
        'Longitude' - Station longitude
        'Elevation_masl' - Station elevation, meters above sea level
        'source' - Source of station data
        'SRid' - Station_event table id
        'event_id' - event id
        'station_id' - Station id (same as Sid)
        'stasource_radius_km' - Distance from event center location to station in km
        'az' - Azimuth from source to station in degrees
        'baz' - Backazimuth (station to source) in degrees
        'detect_HF' - 1 if detected on this station in 1-5 Hz band (VEL), 0 if not, None if unknown
        'detect_LP' - 1 if detected on this station in 20-60 sec band (DISP), 0 if not, None if unknown
        'starttimeHF' - starttime of HF signal on this station in UTCDateTime
        'endtimeHF' - endtime of HF signal on this station in UTCDateTime
        'starttimeHF_kurtosis' - WILL BE DELETED
        'starttimeLP - starttime of LP signal on this station in UTCDateTime
        'endtimeLP' - endtime of LP signal on this station in UTCDateTime
        'peakfreqraw' - Peak frequency of raw signal (highest fourier amplitude)
        'meansqfreqraw' - Mean squared frequency of raw signals
        'meansqfreqSNRraw' - Mean squared frequency of raw signal where SNR greater than 2.0
          relative to before the signal
        'duration_secHF' - Duration of HF signal on this station in seconds
        'duration_secLP'- Duration of LP signal on this station in seconds
        'absmaxampHF' - Maximum absolute value of HF signal in m/s
        'p2pmaxampHF' - Maximum peak to peak value of HF signal in m/s
        'absmaxampLP' - Maximum absolute value of LP signal in m at this station
        'p2pmaxampLP' - Maximum peak to peak value of LP signal in m at this station

    """
    connection = None
    connection = lite.connect(database, detect_types=lite.PARSE_DECLTYPES)
    connection.text_factory = str
    if maxradius is None:
        maxradius = 2000000.
    # Get all station data
    with connection:
        cursor = connection.cursor()
        try:
            event_id = int(event_id)
            cursor_output = (cursor.execute('SELECT * FROM stations s INNER JOIN sta_nearby n ON n.station_id = s.Sid WHERE n.event_id = ? AND n.stasource_radius_km < ? AND n.stasource_radius_km > ? ORDER BY n.stasource_radius_km', (event_id, maxradius, minradius)))
            dat = cursor_output.fetchall()
            temp = cursor.description
            names = [t[0] for t in temp]
            # Turn into list of dictionaries
            staDict = []
            for d in dat:
                staDict.append(dict(list(zip(names, d))))
        except Exception as e:
            print(e)
            return
    # Now filter out undesired entries
    newstas = []
    for stad in staDict:
        if not both:
            if detectHF and (stad['detect_HF'] == 0 or stad['detect_HF'] is None):
                continue
            if detectLP and (stad['detect_LP'] == 0 or stad['detect_LP'] is None):
                continue
            if detectHF is not None:
                if not detectHF and stad['detect_HF'] == 1:
                    continue
            if detectLP is not None:
                if not detectLP and stad['detect_LP'] == 1:
                    continue
        else:  # Get rid of station is both HF and LP are non-detections
            if (stad['detect_HF'] == 0 or stad['detect_HF'] is None) and (stad['detect_LP'] == 0 or stad['detect_LP'] is None):
                continue

        if ampdetectHF is not None:
            if ampdetectHF and stad['absmaxampHF'] is None:
                continue
            if not ampdetectHF and stad['absmaxampHF'] is not None:
                continue
        if ampdetectLP is not None:
            if ampdetectLP and stad['absmaxampLP'] is None:
                continue
            if not ampdetectLP and stad['absmaxampLP'] is not None:
                continue
        if durationHF is not None:
            if durationHF and stad['duration_secHF'] is None:
                continue
            if not durationHF and stad['duration_secHF'] is not None:
                continue
        if durationLP is not None:
            if durationLP and stad['duration_secLP'] is None:
                continue
            if not durationLP and stad['duration_secLP'] is not None:
                continue

        if chanuse == '*':
            pass
        elif stad['Channel'] not in chanuse:
            continue
        try:
            temp = stad['starttimeHF'].split(' ')
            stad['starttimeHF'] = UTCDateTime(temp[0]+'T'+temp[1])
            temp = stad['endtimeHF'].split(' ')
            stad['endtimeHF'] = UTCDateTime(temp[0]+'T'+temp[1])
        except Exception as e:
            pass
            #print(e)
        try:
            temp = stad['starttimeLP'].split(' ')
            stad['starttimeLP'] = UTCDateTime(temp[0]+'T'+temp[1])
            temp = stad['endtimeLP'].split(' ')
            stad['endtimeLP'] = UTCDateTime(temp[0]+'T'+temp[1])
        except Exception as e:
            pass
        # If it made it to here, keep it
        newstas.append(stad)

    if numstas is not None and len(newstas) > 0:
        # Order by distance
        dists = [sta['stasource_radius_km'] for sta in newstas]
        keep = np.sort(np.unique(dists))
        if len(newstas) > numstas:
            keep = keep[0:numstas]
        keepmax = keep.max()
        newstas2 = [sta for sta in newstas if sta['stasource_radius_km'] <= keepmax]
    else:
        newstas2 = newstas

    return newstas2


def attach_distaz(st, event_lat, event_lon, database='/Users/kallstadt/LSseis/landslideDatabase/lsseis.db'):
    """add distance, az, baz for a given event location to each trace in st from database

    Args:
        st: obspy Stream to attach distaz to
        event_lat (float): latitude in decimal degrees of event
        event_lon (float): longitude in decimal degrees of event

    Returns:
        same obspy Stream but with trace.stats.rdist, trace.stats.azimuth, and trace.stats.back_azimuth added for each trace. Events that were not in the database are deleted

    """
    #get station lat lon from database
    connection = None
    connection = lite.connect(database)
    #client = Client()
    stnew = Stream()
    for trace in st:
        chan = trace.stats.channel
        if chan[2] in ('R', 'T'):  # if rotated already, just pretend it's a Z channel
            chan = chan[:2]+'Z'
        with connection:
            cursor = connection.cursor()
            try:
                #get the Sid from stations table
                cursor_output = cursor.execute('SELECT Sid,Latitude,Longitude FROM stations WHERE Name=? AND Channel=?', (trace.stats.station, chan))
            except Exception as e:
                    print(e)
                    continue
            retrieved_data = cursor_output.fetchone()
            if retrieved_data is not None:
                Sid, sta_lat, sta_lon = retrieved_data
            else:
                print((trace.stats.station + '- not in database, skipping'))
                continue
        backazimuth, azimuth, distance = reviewData.pyproj_distaz(sta_lat, sta_lon, event_lat, event_lon)
        #result = client.distaz(sta_lat, sta_lon, event_lat, event_lon)
        trace.stats.rdist = distance  # result['distance']*111.32
        trace.stats.azimuth = azimuth  # result['azimuth']
        trace.stats.back_azimuth = backazimuth  # result['backazimuth']
        stnew += Stream(trace)
    st = stnew
    return st


def attach_coords(st, database='/Users/kallstadt/LSseis/landslideDatabase/lsseis.db'):
    """attach coordinates to stations in Stream st by getting info from IRIS

    Args:
        st: obspy Stream to attach distaz to
        database (str): file path to sqlite3 database to use as source of station coordinates

    Returns:
        same obspy Stream but with coordinates attached.

    """
    from obspy.core import AttribDict
    #get station lat lon from database
    connection = None
    connection = lite.connect(database)
    for trace in st:
        chan = trace.stats.channel
        if chan[2] in ('R', 'T'):  # if rotated already, just pretend it's a Z channel
            chan = chan[:2]+'Z'
        with connection:
            cursor = connection.cursor()
            try:
                #get the Sid from stations table
                cursor_output = cursor.execute('SELECT Latitude, Longitude, Elevation_masl FROM stations WHERE Name=? AND Channel=?', (trace.stats.station, chan))
            except Exception as e:
                    print(e)
                    continue
            retrieved_data = cursor_output.fetchone()
            if retrieved_data is not None:
                lat, lon, elev = retrieved_data
            else:
                print((trace.stats.station + ' not in database, skipping'))
                continue
        trace.stats.coordinates = AttribDict({'latitude': lat, 'longitude': lon, 'elevation': elev})
    return st
