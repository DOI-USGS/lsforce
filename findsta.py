#!/usr/bin/env python

"""
Functions for interacting with and populating the seismogenic landslide database
"""

import sqlite3 as lite
from datetime import datetime
#import urllib2
#from obspy.iris import Client
from obspy import Stream, UTCDateTime
import numpy as np
from reviewData import reviewData


def getEventInfo(event_id, database='/Users/kallstadt/LSseis/landslideDatabase/lsseis.db'):
    """This function pulls the entire entry from the events table into a dictionary (e.g. Lats/Lons, volumes, everything in row)
    :param event_id: id number corresponding to the event of interest in the database
    :type event_id: integer
    :param database: Full file path to SQLite3 landslide database file
    :type database: string
    :returns eventDict: Dictionary containing all entries of all fields from event table
     of database
    :type eventDict: Dictionary

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
            cursor_output = (cursor.execute('SELECT * FROM events WHERE Eid = ?', (event_id, )))
            dat = cursor_output.fetchall()[0]
            temp = cursor.description
            names = [t[0] for t in temp]
            eventDict = dict(zip(names, dat))
            # Convert to the right units
            eventDict['StartTime'] = UTCDateTime(eventDict['StartTime'])
            eventDict['EndTime'] = UTCDateTime(eventDict['EndTime'])
        except Exception as e:
            print(e)
            return

        return eventDict


def getStaInfo(event_id, maxradius=None, minradius=0., detectHF=None, detectLP=None,
               ampdetectHF=None, ampdetectLP=None, durationHF=None, durationLP=None,
               chanuse='*', database='/Users/kallstadt/LSseis/landslideDatabase/lsseis.db'):
    """Get all the station detection information with or without constraints on radius,
    whether or not it was detected in HF or LP bands, and whether or not there
    are amplitude and duration measurements in HF or LP bands.
    :param event_id: id number corresponding to the event of interest in the database
    :type event_id: integer
    :param maxradius: Maximum radius to search in km
    :param minradius: Minimum radius to search in km
    :param detectHF: Return HF detections? None, True, or False
    :param detectLP: Return LP detections? None, True, or False
    :param ampdetectHF: Return stations with HF amplitude measurements? None, True, or False
    :param ampdetectLP: Return stations with HF amplitude measurements? None, True, or False
    :param durationHF: Return stations with HF duration measurements? None, True or False
    :param durationLP: Return stations with LP duration measurements? None, True or False
    :param chanuse: Single string of channel names, return only stations that match channels defined,
     '*' for all channels
    :param database: Full file path to SQLite3 landslide database file
    :type database: string
    :returns eventDict: Dictionary containing all entries of all fields from event table
     of database
    :type eventDict: Dictionary

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
    # Get all station dadta
    with connection:
        cursor = connection.cursor()
        try:
            cursor_output = (cursor.execute('SELECT * FROM stations s INNER JOIN sta_nearby n ON n.station_id = s.Sid WHERE n.event_id = ? AND n.stasource_radius_km < ? AND n.stasource_radius_km > ? ORDER BY n.stasource_radius_km', (event_id, maxradius, minradius)))
            dat = cursor_output.fetchall()
            temp = cursor.description
            names = [t[0] for t in temp]
            # Turn into dictionary using Sid as key
            staDict = {}
            for d in dat:
                staDict[d[0]] = dict(zip(names, d))
        except Exception as e:
            print(e)
            return
    # Now filter out undesired entries
    removekeys = []
    for key in staDict:
        if detectHF is True and (staDict[key]['detect_HF'] == 0 or staDict[key]['detect_HF'] is None):
            removekeys.append(key)
            continue
        if detectHF is False and staDict[key]['detect_HF'] == 1:
            removekeys.append(key)
            continue
        if detectLP is True and (staDict[key]['detect_LP'] == 0 or staDict[key]['detect_LP'] is None):
            removekeys.append(key)
            continue
        if detectLP is False and staDict[key]['detect_LP'] == 1:
            removekeys.append(key)
            continue
        if ampdetectHF is True and staDict[key]['absmaxampHF'] is None:
            removekeys.append(key)
            continue
        if ampdetectHF is False and staDict[key]['absmaxampHF'] is not None:
            removekeys.append(key)
            continue
        if ampdetectLP is True and staDict[key]['absmaxampLP'] is None:
            removekeys.append(key)
            continue
        if ampdetectLP is False and staDict[key]['absmaxampLP'] is not None:
            removekeys.append(key)
            continue
        if durationHF is True and staDict[key]['duration_secHF'] is None:
            removekeys.append(key)
            continue
        if durationHF is False and staDict[key]['duration_secHF'] is not None:
            removekeys.append(key)
            continue
        if durationLP is True and staDict[key]['duration_secLP'] is None:
            removekeys.append(key)
            continue
        if durationLP is False and staDict[key]['duration_secLP'] is not None:
            removekeys.append(key)
            continue
        if chanuse == '*':
            pass
        elif staDict[key]['Channel'] not in chanuse:
            removekeys.append(key)
            continue
        try:
            staDict[key]['starttimeHF'] = UTCDateTime(staDict[key]['starttimeHF'])
            staDict[key]['endtimeHF'] = UTCDateTime(staDict[key]['endtimeHF'])
        except Exception as e:
            print(e)
        try:
            staDict[key]['starttimeLP'] = UTCDateTime(staDict[key]['starttimeLP'])
            staDict[key]['endtimeLP'] = UTCDateTime(staDict[key]['endtimeLP'])
        except Exception as e:
            print(e)

    for k in removekeys:
        staDict.pop(k, None)
    return staDict

def attach_distaz(st, event_lat, event_lon, database='/Users/kallstadt/LSseis/landslideDatabase/lsseis.db'):
    """
    add distance, az, baz for a given event location to each trace in st from database
    USAGE
    st = attach_distaz(st, event_lat, event_lon, database='/Users/kallstadt/LSseis/landslideDatabase/lsseis.db')
    INPUTS
    st - obspy Stream
    event_lat - latitude in decimal degrees of event
    event_lon - longitude in decimal degrees of event
    OUTPUTS
    st - same obspy Stream but with trace.stats.rdist, trace.stats.azimuth, and trace.stats.back_azimuth added for each trace. Events that were not in the database are deleted

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
    """
    attach coordinates to stations in Stream st by getting info from IRIS
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
                print((trace.stats.station + '- not in database, skipping'))
                continue
        trace.stats.coordinates = AttribDict({'latitude': lat, 'longitude': lon, 'elevation': elev})
    return st


# BELOW HERE, ALL FUNCTIONS ARE DEPRECATED, REPLACED WITH TWO FUNCTIONS ABOVE
def get_event_datainfo(event_id, maxradius=50, minradius=0, chanuse='*', database='/Users/kallstadt/LSseis/landslideDatabase/lsseis.db'):
    """
    DEPRECATED FUNCTION
    Get the station and channel names of data corresponding to an event - note that database needs to be populated to maxradius for this event or else this will return nothing or only stations out to maximum radius for which it was populated

    USAGE
    stations, networks, channels, locations, datsource, dist, starttime, endtime, lats, lons = get_event_datainfo(event_id, maxradius=50, minradius=0, channels='*', database='/Users/kallstadt/LSseis/landslideDatabase/lsseis.db')

    INPUTS
    event_id - id of the event in the database
    maxradius - maximum radius in which to search for stations
    minradius - minimum radius in which to search for stations
    chanuse - single string of comma separated channels to include, '*' wildcard finds all channels
    database - full file path to landslide database

    OUTPUTS
    stations - list of station names within specified radius
    networks - list of networks corresponding to station list
    channels - list of channels corresponding to station list
    locations - list of location codes corresponding to station list
    datsource - list of location of data corresponding to station list (e.g. IRIS, NCEDC)
    dist - list of distance in km of each station from event
    starttime - start time of event (datetime, not UTCDateTime format)
    endtime - end time of event (datetime, not UTCDateTime format)
    lats - list of station lats
    lons - list of station lons
    """
    connection = None
    connection = lite.connect(database)
    connection.text_factory = str
    with connection:
        cursor = connection.cursor()
        try:
            cursor_output = (cursor.execute(
                'SELECT Name, Network, Channel, LocationCode, source, stasource_radius_km, Latitude, Longitude FROM sta_nearby, stations ON sta_nearby.station_id = stations.Sid WHERE event_id = ? AND stasource_radius_km < ? AND stasource_radius_km > ? ORDER BY stasource_radius_km', (event_id, maxradius, minradius)))
        except Exception as e:
            print(e)
        retrieved_data = cursor_output.fetchall()
        if chanuse == '*':
            stations = [x[0] for x in retrieved_data]
            networks = [x[1] for x in retrieved_data]
            channels = [x[2] for x in retrieved_data]
            locations = [x[3] for x in retrieved_data]
            datsource = [x[4] for x in retrieved_data]
            dist = [x[5] for x in retrieved_data]
            lats = [x[6] for x in retrieved_data]
            lons = [x[7] for x in retrieved_data]
        else:
            stations = [x[0] for x in retrieved_data if x[2] in chanuse]
            networks = [x[1] for x in retrieved_data if x[2] in chanuse]
            channels = [x[2] for x in retrieved_data if x[2] in chanuse]
            locations = [x[3] for x in retrieved_data if x[2] in chanuse]
            datsource = [x[4] for x in retrieved_data if x[2] in chanuse]
            dist = [x[5] for x in retrieved_data if x[2] in chanuse]
            lats = [x[6] for x in retrieved_data if x[2] in chanuse]
            lons = [x[7] for x in retrieved_data if x[2] in chanuse]
    with connection:
        cursor = connection.cursor()
        try:
            cursor_output2 = (cursor.execute(
                'SELECT Starttime, Endtime FROM events WHERE Eid=?', (event_id,)))
        except Exception as e:
            print(e)
        retrieved_data1 = cursor_output2.fetchall()
        starttime, endtime = retrieved_data1[0]
        try:
            starttime = datetime.strptime(starttime, '%Y-%m-%d %H:%M:%S')
        except Exception:
            starttime = datetime.strptime(starttime, '%Y-%m-%d %H:%M')
        try:
            endtime = datetime.strptime(endtime, '%Y-%m-%d %H:%M:%S')
        except Exception:
            endtime = datetime.strptime(endtime, '%Y-%m-%d %H:%M')

        return stations, networks, channels, locations, datsource, dist, starttime, endtime, lats, lons


def get_event_info(event_id, database='/Users/kallstadt/LSseis/landslideDatabase/lsseis.db'):
    """
    DEPRECATED FUNCTION
    Get basic event info

    USAGE
    event_lat, event_lon, starttime, endtime, datloc = get_event_info(event_id, database='/Users/kallstadt/LSseis/landslideDatabase/lsseis.db')

    INPUTS
    event_id - id of event to retrieve
    database - full file path of location of landslide database

    OUTPUTS
    event_lat - latitude of event (float)
    event_lon - longitude of event (float)
    event_time - start time of event (datetime format, not UTCDateTime)
    datloc - location of data for this event (e.g. IRIS, NCEDC, sac:/filepath) (string)
    """
    connection = None
    connection = lite.connect(database)
    with connection:
        cursor = connection.cursor()
        try:
            cursor_output = (cursor.execute(
                'SELECT Latitude, Longitude, Starttime, Endtime, DatLocation FROM events WHERE Eid=?', (event_id,)))
        except Exception as e:
            print(e)
    retrieved_data = cursor_output.fetchall()
    event_lat, event_lon, starttime, endtime, datloc = retrieved_data[0]
    try:
        starttime = datetime.strptime(starttime, '%Y-%m-%d %H:%M:%S')
    except:
        starttime = datetime.strptime(starttime, '%Y-%m-%d %H:%M')  # return event_time as datetime object
    try:
        endtime = datetime.strptime(endtime, '%Y-%m-%d %H:%M:%S')
    except:
        endtime = datetime.strptime(endtime, '%Y-%m-%d %H:%M')  # return event_time as datetime object
    return event_lat, event_lon, starttime, endtime, datloc


def get_event_info_detail(event_id, database='/Users/kallstadt/LSseis/landslideDatabase/lsseis.db'):
    """
    Get detailed info about event (except station info)

    USAGE
    event_lat, event_lon, locuncert_km, starttime, endtime, name, eventtype, volume, maxdistHF_km, maxdistLP_km, datloc = get_event_info_detail(event_id, database='/Users/kallstadt/LSseis/landslideDatabase/lsseis.db')

    INPUTS
    event_id - id of event to retrieve
    database - full file path of location of landslide database

    OUTPUTS
    event_lat - latitude of event (float)
    event_lon - longitude of event (float)
    locuncert_km - location uncertainty (km)
    starttime - event start time (datetime format, not UTCDateTime)
    endtime - event end time
    name - descriptive name of event (not necessarily unique) (string)
    eventtype - type of event e.g. Rock avalanche (string)
    volume - volume, None if no info, or float
    maxdistHF_km - maximum distance (km) to which high frequencies were observable (0.5 - 5 Hz) (float)
    maxdistLP_km - maximum distance (km) to which long periods were observable (0.05 - 0.2 Hz; 5-20 sec) (float)
    """
    connection = None
    connection = lite.connect(database)
    with connection:
        cursor = connection.cursor()
        try:
            cursor_output = (cursor.execute(
                'SELECT Latitude, Longitude, LocUncert_km, StartTime, EndTime, Name, Type, Volume, maxdistHF_km, maxdistLP_km, DatLocation FROM events WHERE Eid=?', (event_id,)))
        except Exception as e:
            print(e)
    retrieved_data = cursor_output.fetchall()
    event_lat, event_lon, locuncert_km, starttime, endtime, name, eventtype, volume, maxdistHF_km, maxdistLP_km, datloc = retrieved_data[0]
    try:
        starttime = datetime.strptime(starttime, '%Y-%m-%d %H:%M:%S')
    except:
        starttime = datetime.strptime(starttime, '%Y-%m-%d %H:%M')  # return event_time as datetime object
    try:
        endtime = datetime.strptime(endtime, '%Y-%m-%d %H:%M:%S')
    except:
        endtime = datetime.strptime(endtime, '%Y-%m-%d %H:%M')  # return event_time as datetime object
    if volume is not None:
        volume = float(volume)
    if maxdistHF_km is not None:
        try:
            maxdistHF_km = float(maxdistHF_km)
        except:
            pass
    if maxdistLP_km is not None:
        try:
            maxdistLP_km = float(maxdistLP_km)
        except:
            pass
    return event_lat, event_lon, locuncert_km, starttime, endtime, name, eventtype, volume, maxdistHF_km, maxdistLP_km, datloc


def get_sta_detect(event_id, detectHF=True, detectLP=None, mindistkm=0., maxdistkm=None, database='/Users/kallstadt/LSseis/landslideDatabase/lsseis.db'):
    """
    Get list of stations that registered an event either in HF or LP domain or both within a given radius
    detectHF and detectLP can be True, False or None (for NULL)
    """
    if detectHF is True:
        detectHF = 1
    elif detectHF is False:
        detectHF = 0
    if detectLP is True:
        detectLP = 1
    elif detectLP is False:
        detectLP = 0
    if maxdistkm is None:
        maxdistkm = 10000000.
    # Connect to database and query it
    connection = None
    connection = lite.connect(database)
    connection.text_factory = str
    with connection:
        cursor = connection.cursor()
        # Find stations that meet criteria
        try:
            if detectLP is None:
                cursor_output = (cursor.execute(
                                 'SELECT s.Name, s.Network, s.Channel, s.LocationCode, s.source, n.SRid, s.Latitude, s.Longitude, n.stasource_radius_km, n.az, n.starttimeHF, n.endtimeHF FROM stations s INNER JOIN sta_nearby n ON n.station_id = s.Sid WHERE n.event_id = ? AND n.stasource_radius_km < ? AND n.stasource_radius_km > ? AND n.detect_HF = ? ORDER BY n.stasource_radius_km', (event_id, maxdistkm, mindistkm, detectHF)))
            elif detectHF is None:
                cursor_output = (cursor.execute(
                                 'SELECT s.Name, s.Network, s.Channel, s.LocationCode, s.source, n.SRid, s.Latitude, s.Longitude, n.stasource_radius_km, n.az, n.starttimeHF, n.endtimeHF FROM stations s INNER JOIN sta_nearby n ON n.station_id = s.Sid WHERE n.event_id = ? AND n.stasource_radius_km < ? AND n.stasource_radius_km > ? AND n.detect_LP = ? ORDER BY n.stasource_radius_km', (event_id, maxdistkm, mindistkm, detectLP)))
            else:
                cursor_output = (cursor.execute(
                                 'SELECT s.Name, s.Network, s.Channel, s.LocationCode, s.source, n.SRid, s.Latitude, s.Longitude, n.stasource_radius_km, n.az, n.starttimeHF, n.endtimeHF FROM stations s INNER JOIN sta_nearby n ON n.station_id = s.Sid WHERE n.event_id = ? AND n.stasource_radius_km < ? AND n.stasource_radius_km > ? AND n.detect_HF = ? AND n.detect_LP = ? ORDER BY n.stasource_radius_km', (event_id, maxdistkm, mindistkm, detectHF, detectLP)))
        except Exception as e:
            print(e)
        retrieved_data = cursor_output.fetchall()

        stations = [x[0] for x in retrieved_data]
        networks = [x[1] for x in retrieved_data]
        channels = [x[2] for x in retrieved_data]
        locations = [x[3] for x in retrieved_data]
        datsource = [x[4] for x in retrieved_data]
        SRid = [x[5] for x in retrieved_data]
        latitudes = [x[6] for x in retrieved_data]
        longitudes = [x[7] for x in retrieved_data]
        dist = [x[8] for x in retrieved_data]
        az = [x[9] for x in retrieved_data]
        stt = [x[10] for x in retrieved_data]
        endt = [x[11] for x in retrieved_data]

        starttimesHF = []
        endtimesHF = []

        for i in np.arange(len(stt)):
            try:
                try:
                    starttimesHF.append(datetime.strptime(stt[i], '%Y-%m-%d %H:%M:%S'))
                except:
                    starttimesHF.append(datetime.strptime(stt[i], '%Y-%m-%d %H:%M'))  # return event_time as datetime object
                try:
                    endtimesHF.append(datetime.strptime(endt[i], '%Y-%m-%d %H:%M:%S'))
                except:
                    endtimesHF.append(datetime.strptime(endt[i], '%Y-%m-%d %H:%M'))  # return event_time as datetime object
            except:
                endtimesHF.append(None)
                starttimesHF.append(None)

        for i, loc in enumerate(locations):
            if '--' in loc:
                locations[i] = ''

    return stations, networks, channels, locations, datsource, SRid, latitudes, longitudes, dist, az, starttimesHF, endtimesHF


def get_ampdetect(event_id, mindistkm=0., maxdistkm=None, database='/Users/kallstadt/LSseis/landslideDatabase/lsseis.db'):
    """
    Get list of stations (and their lats and lons) that have amplitude measurements for a given event
    """
    if maxdistkm is None:
        maxdistkm = 10000000.
    # Connect to database and query it
    connection = None
    connection = lite.connect(database)
    connection.text_factory = str
    with connection:
        cursor = connection.cursor()
        # Find stations that meet criteria
        try:
            cursor_output = (cursor.execute(
                'SELECT s.Name, s.Network, s.Channel, s.LocationCode, s.Latitude, s.Longitude, n.stasource_radius_km, n.az, n.absmaxampHF, n.duration_secHF FROM stations s INNER JOIN sta_nearby n ON n.station_id = s.Sid WHERE n.event_id = ? AND n.stasource_radius_km < ? AND n.stasource_radius_km > ? AND n.absmaxampHF is not NULL ORDER BY n.stasource_radius_km', (event_id, maxdistkm, mindistkm)))
        except Exception as e:
            print(e)
        retrieved_data = cursor_output.fetchall()

        stations = [x[0] for x in retrieved_data]
        networks = [x[1] for x in retrieved_data]
        channels = [x[2] for x in retrieved_data]
        locations = [x[3] for x in retrieved_data]
        latitudes = [x[4] for x in retrieved_data]
        longitudes = [x[5] for x in retrieved_data]
        dist = [x[6] for x in retrieved_data]
        az = [x[7] for x in retrieved_data]
        absamps = [x[8] for x in retrieved_data]
        duration = [x[9] for x in retrieved_data]

    return stations, networks, channels, locations, latitudes, longitudes, dist, az, absamps, duration


def get_ampdetectLP(event_id, mindistkm=0., maxdistkm=None, database='/Users/kallstadt/LSseis/landslideDatabase/lsseis.db'):
    """
    Get list of stations (and their lats and lons) that have amplitude measurements for a given event
    """
    if maxdistkm is None:
        maxdistkm = 10000000.
    # Connect to database and query it
    connection = None
    connection = lite.connect(database)
    connection.text_factory = str
    with connection:
        cursor = connection.cursor()
        # Find stations that meet criteria
        try:
            cursor_output = (cursor.execute(
                'SELECT s.Name, s.Network, s.Channel, s.LocationCode, s.Latitude, s.Longitude, n.stasource_radius_km, n.az, n.p2pmaxampLP, n.duration_secLP FROM stations s INNER JOIN sta_nearby n ON n.station_id = s.Sid WHERE n.event_id = ? AND n.stasource_radius_km < ? AND n.stasource_radius_km > ? AND n.p2pmaxampLP is not NULL ORDER BY n.stasource_radius_km', (event_id, maxdistkm, mindistkm)))
        except Exception as e:
            print(e)
        retrieved_data = cursor_output.fetchall()

        stations = [x[0] for x in retrieved_data]
        networks = [x[1] for x in retrieved_data]
        channels = [x[2] for x in retrieved_data]
        locations = [x[3] for x in retrieved_data]
        latitudes = [x[4] for x in retrieved_data]
        longitudes = [x[5] for x in retrieved_data]
        dist = [x[6] for x in retrieved_data]
        az = [x[7] for x in retrieved_data]
        absamps = [x[8] for x in retrieved_data]
        duration = [x[9] for x in retrieved_data]

    return stations, networks, channels, locations, latitudes, longitudes, dist, az, absamps, duration


def get_durdetect(event_id, mindistkm=0., maxdistkm=None, database='/Users/kallstadt/LSseis/landslideDatabase/lsseis.db'):
    """
    Get list of stations (and their lats and lons) that have duration measurements for a given event
    """
    if maxdistkm is None:
        maxdistkm = 10000000.
    # Connect to database and query it
    connection = None
    connection = lite.connect(database)
    connection.text_factory = str
    with connection:
        cursor = connection.cursor()
        # Find stations that meet criteria
        try:
            cursor_output = (cursor.execute(
                'SELECT s.Name, s.Network, s.Channel, s.LocationCode, s.Latitude, s.Longitude, n.stasource_radius_km, n.az, n.absmaxampHF, n.duration_secHF, n.meansqfreqSNRraw FROM stations s INNER JOIN sta_nearby n ON n.station_id = s.Sid WHERE n.event_id = ? AND n.stasource_radius_km < ? AND n.stasource_radius_km > ? AND n.duration_secHF is not NULL ORDER BY n.stasource_radius_km', (event_id, maxdistkm, mindistkm)))
        except Exception as e:
            print(e)
        retrieved_data = cursor_output.fetchall()

        stations = [x[0] for x in retrieved_data]
        networks = [x[1] for x in retrieved_data]
        channels = [x[2] for x in retrieved_data]
        locations = [x[3] for x in retrieved_data]
        latitudes = [x[4] for x in retrieved_data]
        longitudes = [x[5] for x in retrieved_data]
        dist = [x[6] for x in retrieved_data]
        az = [x[7] for x in retrieved_data]
        absamps = [x[8] for x in retrieved_data]
        duration = [x[9] for x in retrieved_data]
        meansqfreqSNR = [x[10] for x in retrieved_data]

    return stations, networks, channels, locations, latitudes, longitudes, dist, az, absamps, duration, meansqfreqSNR
