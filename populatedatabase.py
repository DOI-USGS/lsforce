#!/usr/bin/env python

"""
Functions for populating and updating the database
"""

import findsta
from reviewData import reviewData
import numpy as np
import sqlite3 as lite


def initial_populate(event_ids, maxradius=500., IRIS=True, NCEDC=False):
    """
    :param event_ids: either single integer or list or numpy array of integers, eg np.arange(62,65)
    :param maxradius: Maximum radius to go out to in searching for nearby stations in stasource_radius_km
    :param IRIS: If True, will search IRIS for stations nearby
    :param NCEDC: If True, will search NCEDC for stations nearby
    """
    for event_id in event_ids:
        if IRIS is True:
            evdict = findsta.getEventInfo(event_id)
            lines, source = reviewData.get_stations_iris(evdict['Latitude'], evdict['Longitude'],
                                                         evdict['StartTime'], maxradiuskm=maxradius)
            findsta.populate_station_tables(lines, source)
            findsta.populate_station_event_table(event_id, lines, update=True)
            lines = None
            source = None
        if NCEDC is True:
            lines, source = reviewData.get_stations_ncedc(evdict['Latitude'], evdict['Longitude'],
                                                          evdict['StartTime'], maxradiuskm=maxradius)
            findsta.populate_station_tables(lines, source)
            findsta.populate_station_event_table(event_id, lines, update=True)
            lines = None
            source = None


def recalculate_distances(event_ids, database='/Users/kallstadt/LSseis/landslideDatabase/lsseis.db'):
    """Recalculate source to station distances. Should be done if lat/lons were updated for either
    source or stations in the table
    :param event_ids: either single integer or list or numpy array of integers, eg np.arange(62,65)
    """
    connection = None
    connection = lite.connect(database)
    with connection:
        for event_id in event_ids:
        # Get all the SRids and corresponding station_id's for event_id
            evdict = findsta.getEventInfo(event_id)
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
                backazimuth, azimuth, distance = reviewData.pyproj_distaz(sta_lat, sta_lon, evdict['Latitude'], evdict['Longitude'])
                cursor.execute("UPDATE sta_nearby SET stasource_radius_km = ?, az = ?, baz = ? WHERE SRid=?", ('%3.3f' % distance, '%3.2f' % azimuth, '%3.2f' % backazimuth, SRid[i]))


def populate_station_tables(lines, source, database='/Users/kallstadt/LSseis/landslideDatabase/lsseis.db'):
    """
    Helper code used to put the station info into the stations table of database if it's not already there
    """
    #loop over each line, parse info, see if it already is in table, if not, insert it
    connection = None
    connection = lite.connect(database)
    with connection:
        cursor = connection.cursor()
        val = 0
        val2 = 0
        for line in lines:
            try:
                #see if the station is already there
                cursor_output = (cursor.execute(
                    'SELECT * FROM stations WHERE Name=? AND Channel=? AND Network=?', (line[1], line[3], line[0])))
            except Exception as e:
                print(e)
            retrieved_data = cursor_output.fetchall()
            if len(retrieved_data) == 0:  # if the station isn't already there, put it there
                try:
                    cursor.execute(
                        'INSERT INTO stations(Network,Name,LocationCode,Channel,Latitude,Longitude,Elevation_masl,source) VALUES(?,?,?,?,?,?,?,?)',
                        (line[0], line[1], line[2], line[3], line[4], line[5], line[6], source))
                    val += 1
                except Exception as f:
                    print(f)
            else:
                val2 += 1
        print(('added %s entries to stations table, %s already were there' % (val, val2)))


def populate_station_event_table(event_id, lines, database='/Users/kallstadt/LSseis/landslideDatabase/lsseis.db', update=True):
    """
    Also, put the station in the sta_nearby table for this event and calculate distance, az, baz
    """
    #get event lat lon from database
    connection = None
    connection = lite.connect(database)
    with connection:
        cursor = connection.cursor()
        try:
            cursor_output = cursor.execute('SELECT Latitude, Longitude FROM events WHERE Eid=?', (event_id,))
        except Exception as e:
            print(e)
            return
    retrieved_data = cursor_output.fetchall()
    event_lat, event_lon = retrieved_data[0]

    #get Sid and lat lon for each entry in line
    val = 0
    val1 = 0
    val2 = 0
    valbad = 0
    for line in lines:
        with connection:
            try:
                #get the Sid from stations table
                cursor_output = cursor.execute(
                    'SELECT Sid,Latitude,Longitude FROM stations WHERE Name=? AND Channel=? AND Network=?', (line[1], line[3], line[0]))
                retrieved_data = cursor_output.fetchone()
                Sid, sta_lat, sta_lon = retrieved_data
            except:
                sta_lat = None
                sta_lon = None
            try:
                cursor_output = cursor.execute('SELECT SRid FROM sta_nearby WHERE event_id=? AND station_id=?', (event_id, Sid))
                retrieved_data = cursor_output.fetchall()
            except:
                retrieved_data = []
        if len(retrieved_data) == 0:
            #if it isn't already there, use iris's distaz webservice to calculate stuff and save it in the table
            #result = client.distaz(sta_lat, sta_lon, event_lat, event_lon)
            backazimuth, azimuth, distance = reviewData.pyproj_distaz(sta_lat, sta_lon, event_lat, event_lon)
            with connection:
                try:
                    #cursor.execute('INSERT INTO sta_nearby(event_id,station_id,stasource_radius_km,az,baz) VALUES(?,?,?,?,?)', (event_id, Sid, result['distance']*111.32, result['azimuth'], result['backazimuth']))
                    cursor.execute('INSERT INTO sta_nearby(event_id,station_id,stasource_radius_km,az,baz) VALUES(?,?,?,?,?)', (event_id, Sid, '%3.3f' % distance, '%3.2f' % azimuth, '%3.2f' % backazimuth))
                    val += 1
                except Exception as f:
                    print(f)
                    valbad += 1
        elif update is True:
            #if it isn't already there, use iris's distaz webservice to calculate stuff and save it in the table
            with connection:
                try:
                    #result = client.distaz(sta_lat, sta_lon, event_lat, event_lon)
                    backazimuth, azimuth, distance = reviewData.pyproj_distaz(sta_lat, sta_lon, event_lat, event_lon)
                    #cursor.execute("UPDATE sta_nearby SET stasource_radius_km = ?, az = ?, baz = ? WHERE event_id =? AND station_id =?;",(result['distance']*111.32, result['azimuth'], result['backazimuth'], event_id, Sid))
                    cursor.execute("UPDATE sta_nearby SET stasource_radius_km = ?, az = ?, baz = ? WHERE event_id =? AND station_id =?;", ('%3.3f' % distance, '%3.2f' % azimuth, '%3.2f' % backazimuth, event_id, Sid))
                    val1 += 1
                except Exception as f:
                    print(f)
                    valbad += 1
        else:
            val2 += 1
    print(('added %s entries to sta_nearby table, updated %s entries, %s left as is, %s not added because of error' % (val, val1, val2, valbad)))
