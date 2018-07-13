#!/usr/bin/env python

"""
Functions for cleaning up the database, renaming files and so on in preparation for data release
"""
import numpy as np
import sqlite3 as lite
from shutil import copy, rmtree
import datetime
import os
import csv
import findsta
from PIL import Image
from PIL.ExifTags import TAGS, GPSTAGS
import glob
import time
import reverse_geocoder as rg


def get_exif(imgfilename):
    """Returns info from exif data of an Image, if it exists.
    Modified from https://gist.github.com/erans/983821#file-get_lat_lon_exif_pil-py-L4

    Args:
        imgfilename (str): path to image filenames

    Returns:
        tuple: tuple (latitude, longitude, date) where,
            * latitude: latitude of where photo was taken from, in degrees, from GPSTAGS
            * longitude: longitude of where photo was taken from, in degrees, from GPSTAGS
            * date: formatted string of date and time image was taken

    """
    try:
        image = Image.open(imgfilename)
        exif_data = {}
        info = image._getexif()
        if info:
            for tag, value in list(info.items()):
                decoded = TAGS.get(tag, tag)
                if decoded == "GPSInfo":
                    gps_data = {}
                    for t in value:
                        sub_decoded = GPSTAGS.get(t, t)
                        gps_data[sub_decoded] = value[t]

                    exif_data[decoded] = gps_data
                else:
                    exif_data[decoded] = value

            if 'DateTime' in exif_data:
                # convert to datetime object
                pass
            else:
                date = None

            if 'GPSInfo' in exif_data:
                lat = exif_data['GPSInfo']['GPSLatitude']
                latitude = lat[0][0]/lat[0][1] + (lat[1][0]/lat[1][1])/60. + (lat[2][0]/lat[2][1])/3600.
                lon = exif_data['GPSInfo']['GPSLongitude']
                longitude = lon[0][0]/lon[0][1] + (lon[1][0]/lon[1][1])/60. + (lon[2][0]/lon[2][1])/3600.
                if exif_data['GPSInfo']['GPSLatitudeRef'] is 'S':
                    latitude = -latitude
                if exif_data['GPSInfo']['GPSLongitudeRef'] is 'W':
                    longitude = -longitude
                #dat = exif_data['GPSInfo']['GPSDateStamp'].split(':')
                #tim = exif_data['GPSInfo']['GPSTimeStamp']
                #date = datetime.datetime(dat[0], dat[1], dat[2], tim[0][0]/tim[0][1], tim[1][0]/tim[1][1], tim[2][0]/tim[2][1])
                #date = date.strftime('%Y-%m-%d %H:%M:%S')
            else:
                latitude = None
                longitude = None
            if 'DateTime' in exif_data:
                # convert to datetime object
                dat, tim = exif_data['DateTime'].split(' ')
                dat = dat.split(':')
                dat = [int(d) for d in dat]
                tim = tim.split(':')
                date = datetime.datetime(int(dat[0]), int(dat[1]), int(dat[2]), int(tim[0]), int(tim[1]), int(tim[2]))
                date = date.strftime('%Y-%m-%d %H:%M:%S')
            else:
                date = None

        else:
            #print('No exif data returned')
            latitude = None
            longitude = None
            date = None
    except Exception:  # as e:
        #print('No exif data returned: %s' % e)
        latitude = None
        longitude = None
        date = None

    return latitude, longitude, date


def prepare_new(newdbname=None, newifname=None, database='/Users/kallstadt/LSseis/landslideDatabase/lsseis.db',
                infofiles='/Users/kallstadt/LSseis/landslideDatabase/InfoFiles', shortfilen=None):
    """Setup to make a cleaned copy of the database along with flat files for each event

    Args:
        newdbname (str): path and file name for new copy of database, if None, will use current name (database)
            with date appended
        newifname (dict): directory where new copies of info files should be placed, if None, will use current name
            (infofiles) with date appended
        database (str): path to current database file
        infofiles (str): path to directory where current info files are located
        shortfilen (str): relative file path to append to new file locations when they are replaced in database, optional

    Returns:
        tuple (eids, relpath, newdbname, newifname, shortfilen) where,
            * eids: list of event ids
            * relpath: directory where infofiles is located
            * newdbname: name of new database
            * newifname: name of new directory for info files
            * shortfilen: relative file path to append to new file locations when they are replaced in database

    """
    relpath = os.path.dirname(infofiles)
    # Create new copy of database and new copy of InfoFiles
    if newdbname is None:
        time1 = datetime.datetime.utcnow().strftime('%d%b%Y')
        fn, ext = os.path.splitext(database)
        # add datetime to end
        newdbname = '%s_%s%s' % (fn, time1, ext.lower())
    try:
        os.remove(newdbname)
    except Exception as e:
        print(e)

    try:
        rmtree(os.path.join(os.path.dirname(newdbname), 'references'))
    except Exception as e:
        print(e)

    try:
        os.remove(os.path.join(os.path.dirname(newdbname), 'Events.csv'))
    except Exception as e:
        print(e)

    try:
        os.remove(os.path.join(os.path.dirname(newdbname), 'references.csv'))
    except Exception as e:
        print(e)

    if newifname is None:
        time1 = datetime.datetime.utcnow().strftime('%d%b%Y')
        fn, ext = os.path.splitext(database)
        # add datetime to end
        newifname = '%s_%s' % (infofiles, time1)
    try:
        rmtree(newifname)
    except Exception as e:
        print(e)

    copy(database, newdbname)

    try:
        os.mkdir(newifname)
    except Exception as e1:
        print(e1)

    # Make top level reference folder
    refdir = os.path.join(os.path.dirname(newdbname), 'references')
    try:
        os.mkdir(refdir)
    except Exception as e:
        print(e)

    if shortfilen is None:
        shortfilen = newifname.split('landslideDatabase/')[1]

    eids = get_eids(newdbname)

    for e1 in eids:
        # Make new InfoFiles folders for each event
        evinf = findsta.getEventInfo(e1, database=newdbname)
        evname = evinf['Name'].strip().replace(' ', '').replace(',', '')
        evdat = evinf['StartTime'].strftime('%d%b%Y')
        fulldir = os.path.join(newifname, '%s_%s_%s' % (e1, evname, evdat))
        try:
            os.mkdir(fulldir)
        except Exception as e:
            print(e)

    return eids, relpath, newdbname, newifname, shortfilen


def get_eids(newdbname):
    """Get all event ids in database

    Args:
        newdbname (str): full file extension to database to read

    Returns:
        list of event_ids

    """
    # Connect to database
    connection = None
    connection = lite.connect(newdbname)

    # Get list of event ids
    with connection:
        cursor = connection.cursor()
        cursor_output = cursor.execute("""SELECT Eid FROM events""")
    eid1 = cursor_output.fetchall()
    eids = [e2[0] for e2 in eid1]
    eids.reverse()
    connection.close()
    return eids


def clean_stanearby(eids, database=None):
    """Set all stations beyond max detection limits to NULL and detect_LP to NULL for any short period stations

    Args:
        eids (list): list of event_ids
        database (str): full file extension to sqlite3 database to clean

    """
    for e1 in eids:
        evinf = findsta.getEventInfo(e1, database=database)
        # Connect to database
        connection = None
        connection = lite.connect(database)
        cursor = connection.cursor()

        # pull all info from sta_nearby
        cursor_output = cursor.execute("""SELECT SRid, stasource_radius_km FROM sta_nearby WHERE event_id=?""", (e1,))
        stadata = cursor_output.fetchall()
        for temp in stadata:
            SRid, stasource_radius_km = temp
            if stasource_radius_km > evinf['maxdistHF_km']:
                with connection:
                    connection.execute('UPDATE sta_nearby SET detect_HF=? WHERE SRid = ?', (None, SRid))
            if stasource_radius_km > evinf['maxdistLP_km']:
                with connection:
                    connection.execute('UPDATE sta_nearby SET detect_LP=? WHERE SRid = ?', (None, SRid))
            # else:
            #     pt = cursor.execute("""SELECT Channel FROM sta_nearby, stations ON sta_nearby.station_id = stations.Sid WHERE sta_nearby.SRid=?""", (SRid,))
            #     channel = pt.fetchone()
            #     if channel[:2] == 'EH':
            #         with connection:
            #             connection.execute('UPDATE sta_nearby SET detect_LP=? WHERE SRid = ?', (None, SRid))
        connection.close()
        time.sleep(0.1)


def clean_photos(eids, relpath, newdbname, newifname, shortfilen):
    """Make clean copies of all photos/figures and give them uniform names, modify database to reflect changes,
    create flat file for each event summarizing photos.

    Args:
        eids (list): list of event ids
        relpath (str): directory where infofiles is located
        newdbname (str): name of new database
        newifname (str): name of new directory for info files
        shortfilen (str): relative file path to append to new file locations when they are replaced in database
    """

    # Go through each event
    for e1 in eids:
        # Get basic event information
        evinf = findsta.getEventInfo(e1, database=newdbname)
        evname = evinf['Name'].strip().replace(' ', '').replace(',', '')
        evdat = evinf['StartTime'].strftime('%d%b%Y')
        fulldir = os.path.join(newifname, '%s_%s_%s' % (e1, evname, evdat))

        # Connect to database
        connection = None
        connection = lite.connect(newdbname)
        cursor = connection.cursor()

        # Get max photo id currently existing
        idlist = cursor.execute("""SELECT photo_id FROM photos""")
        phoids = idlist.fetchall()
        phoids = [id1[0] for id1 in phoids]
        currentpid = np.max(phoids) + 1

        # Go through and extract and rename all of the associated files, make summary file on each page
        cursor_output = cursor.execute("""SELECT photo_id FROM photos WHERE event_id=?""", (e1,))
        if cursor_output.fetchone() is not None:
            cursor_output = cursor.execute("""SELECT photo_id, photographer, description, file_extension, apply_also FROM photos WHERE event_id=?""", (e1,))
            photodata = cursor_output.fetchall()
            photodir = os.path.join(fulldir, 'photos_figures')
            try:
                os.mkdir(photodir)
            except Exception as e:
                print(e)

            csvfilename = os.path.join(fulldir, '%s_%s_%s_photos_figures.csv' % (e1, evname, evdat))
            if os.path.isfile(csvfilename):
                perm = 'a'
            else:
                perm = 'wb'
            with open(csvfilename, perm) as csvfile:
                num = 1
                # Put something in for the first line
                writer = csv.writer(csvfile)
                if perm == 'wb':
                    writer.writerow(['Photos of %s %s %s (Event_id %d)' % (evinf['StartTime'].strftime('%d%b%YT%H:%M'), evinf['Name'], evinf['Type'], evinf['Eid'])])
                    writer.writerow(['Note: Some of these files may be located in the folders for other events that share background information'])
                    writer.writerow(['photo_id', 'photographer', 'description', 'date', 'latitude', 'longitude', 'file_extension'])
                for temp in photodata:
                    # go through each entry, extract names, rename, copy over, write it down in csv
                    photo_id, photographer, description, file_extension, apply_also = temp
                    photoshort = photographer.replace(' ', '').replace(',', '-').strip('.')
                    # get list of files if there is more than one
                    if ';' in file_extension:
                        manyfiles = file_extension.replace('; ', ';').split(';')
                        filenames = [os.path.join(relpath, mf) for mf in manyfiles]
                        flag = 1
                        fullpath = 'none'
                    else:
                        flag = 0
                        fullpath = os.path.join(relpath, file_extension)

                    if os.path.isdir(fullpath) or flag == 1:
                        if flag == 0:
                            filenames = glob.glob(os.path.join(fullpath, '*'))
                        cnt = 1
                        for filen in filenames:
                            latitude, longitude, date = get_exif(filen)
                            fn, ext = os.path.splitext(filen)
                            basename = os.path.basename(filen)
                            if '_waveforms_' in basename:
                                name, ext1 = os.path.splitext(basename)
                                name = name.split('_waveforms')
                                name = 'waveforms' + name[-1]
                                temp = '%s_%s_%s_%s' % (e1, evname, evdat, name)
                            elif '_spectr' in basename:
                                name, ext1 = os.path.splitext(basename)
                                name = name.split('_spectr')
                                name = 'spectr' + name[-1]
                                temp = '%s_%s_%s_%s' % (e1, evname, evdat, name)
                            elif '_displspectr' in basename:
                                name, ext1 = os.path.splitext(basename)
                                name = name.split('_displspectr')
                                name = 'displspectr' + name[-1]
                                temp = '%s_%s_%s_%s' % (e1, evname, evdat, name)
                            elif '_teles' in basename:
                                name, ext1 = os.path.splitext(basename)
                                name = name.split('_teles')
                                name = '_teles' + name[-1]
                                temp = '%s_%s_%s_%s' % (e1, evname, evdat, name)
                            else:
                                temp = '%s_%s_%s_source-%s_%d' % (e1, evname, evdat, photoshort, num)
                                num += 1
                            newfilen = os.path.join(fulldir, 'photos_figures', temp + ext.lower())
                            newfilenrel = os.path.join(shortfilen, '%s_%s_%s' % (e1, evname, evdat), 'photos_figures', temp + ext.lower())
                            copy(filen, newfilen)
                            if cnt == 1:
                                # replace original entry
                                # try:
                                with connection:
                                    connection.execute('UPDATE photos SET latitude=?, longitude=?, date = ?,  file_extension=? WHERE photo_id = ?', (latitude, longitude, date, newfilenrel, photo_id))
                                writer.writerow([photo_id, photographer, description, date, latitude, longitude, newfilenrel])
                                # except Exception as e:
                                #     print e

                            else:
                                # add new entry
                                # try:
                                with connection:
                                    connection.execute('INSERT INTO photos(photo_id, event_id, photographer, description, latitude, longitude, date, file_extension, apply_also) VALUES (?,?,?,?,?,?,?,?,?)',
                                                       (currentpid, e1, photographer, description, latitude, longitude, date, newfilenrel, apply_also))
                                writer.writerow([currentpid, photographer, description, date, latitude, longitude, newfilenrel])
                                # except Exception as e:
                                #     print e
                                currentpid += 1
                            cnt += 1

                            if apply_also is not None:
                                aids = apply_also.replace(' ', '').split(',')
                                aids = [int(a1) for a1 in aids]
                                for a1 in aids:
                                    # Try to make file, if not already there
                                    aevinf = findsta.getEventInfo(a1, database=newdbname)
                                    aevname = aevinf['Name'].strip().replace(' ', '').replace(',', '')
                                    aevdat = aevinf['StartTime'].strftime('%d%b%Y')
                                    afulldir = os.path.join(newifname, '%s_%s_%s' % (a1, aevname, aevdat))
                                    acsvfilename = os.path.join(afulldir, '%s_%s_%s_photos_figures.csv' % (a1, aevname, aevdat))
                                    if os.path.isfile(acsvfilename):
                                        perm = 'a'
                                    else:
                                        perm = 'wb'
                                    with open(acsvfilename, perm) as acsvfile:
                                        awriter = csv.writer(acsvfile)
                                        if perm == 'wb':
                                            # if new Put something in for the first line
                                            awriter.writerow(['Photos of %s %s %s (Event_id %d)' % (aevinf['StartTime'].strftime('%d%b%YT%H:%M'), aevinf['Name'], aevinf['Type'], aevinf['Eid'])])
                                            awriter.writerow(['Note: Some of these files may be located in the folders for other events that share background information'])
                                            awriter.writerow(['photo_id', 'photographer', 'description', 'date', 'latitude', 'longitude', 'file_extension'])
                                        awriter.writerow([currentpid, photographer, description, date, latitude, longitude, newfilenrel])

                    else:
                        latitude, longitude, date = get_exif(fullpath)
                        fn, ext = os.path.splitext(fullpath)
                        basename = os.path.basename(fullpath)
                        if '_waveforms_' in basename:
                            name, ext1 = os.path.splitext(basename)
                            name = name.split('_waveforms')
                            name = 'waveforms' + name[-1]
                            temp = '%s_%s_%s_%s' % (e1, evname, evdat, name)
                        elif '_spectr' in basename:
                            name, ext1 = os.path.splitext(basename)
                            name = name.split('_spectr')
                            name = 'spectr' + name[-1]
                            temp = '%s_%s_%s_%s' % (e1, evname, evdat, name)
                        elif '_displspectr' in basename:
                            name, ext1 = os.path.splitext(basename)
                            name = name.split('_displspectr')
                            name = 'displspectr' + name[-1]
                            temp = '%s_%s_%s_%s' % (e1, evname, evdat, name)
                        elif '_teles' in basename:
                            name, ext1 = os.path.splitext(basename)
                            name = name.split('_teles')
                            name = '_teles' + name[-1]
                            temp = '%s_%s_%s_%s' % (e1, evname, evdat, name)
                        else:
                            temp = '%s_%s_%s_source-%s_%d' % (e1, evname, evdat, photoshort, num)
                            num += 1
                        newfilen = os.path.join(fulldir, 'photos_figures', temp + ext.lower())
                        newfilenrel = os.path.join(shortfilen, '%s_%s_%s' % (e1, evname, evdat), 'photos_figures', temp + ext.lower())
                        copy(fullpath, newfilen)
                        num += 1
                        # replace original entry
                        # try:
                        with connection:
                            connection.execute('UPDATE photos SET latitude=?, longitude=?, date = ?,  file_extension=? WHERE photo_id = ?', (latitude, longitude, date, newfilenrel, photo_id))
                            writer.writerow([photo_id, photographer, description, date, latitude, longitude, newfilenrel])
                        # except Exception as e:
                        #     print e
                        if apply_also is not None:
                            aids = apply_also.replace(' ', '').split(',')
                            aids = [int(a1) for a1 in aids]
                            for a1 in aids:
                                # Try to make file, if not already there
                                aevinf = findsta.getEventInfo(a1, database=newdbname)
                                aevname = aevinf['Name'].strip().replace(' ', '').replace(',', '')
                                aevdat = aevinf['StartTime'].strftime('%d%b%Y')
                                afulldir = os.path.join(newifname, '%s_%s_%s' % (a1, aevname, aevdat))
                                acsvfilename = os.path.join(afulldir, '%s_%s_%s_photos_figures.csv' % (a1, aevname, aevdat))
                                if os.path.isfile(acsvfilename):
                                    perm = 'a'
                                else:
                                    perm = 'wb'
                                with open(acsvfilename, perm) as acsvfile:
                                    awriter = csv.writer(acsvfile)
                                    if perm == 'wb':
                                        # if new Put something in for the first line
                                        awriter.writerow(['Photos of %s %s %s (Event_id %d)' % (aevinf['StartTime'].strftime('%d%b%YT%H:%M'), aevinf['Name'], aevinf['Type'], aevinf['Eid'])])
                                        awriter.writerow(['Note: Some of these files may be located in the folders for other events that share background information'])
                                        awriter.writerow(['photo_id', 'photographer', 'description', 'date', 'latitude', 'longitude', 'file_extension'])
                                    awriter.writerow([currentpid, photographer, description, date, latitude, longitude, newfilenrel])

    connection.close()
    time.sleep(2.0)


def clean_gis(eids, relpath, newdbname, newifname, shortfilen):
    """Make clean copies of all gis/map figures and give them uniform names, modify database to reflect changes,
        create flat file for each event summarizing the gis files.

    Args:
        eids (list): list of event ids
        relpath (str): directory where infofiles is located
        newdbname (str): name of new database
        newifname (str): name of new directory for info files
        shortfilen (str): relative file path to append to new file locations when they are replaced in database

    """
    # Go through each event
    for e1 in eids:
        # Get basic event information
        evinf = findsta.getEventInfo(e1, database=newdbname)
        evname = evinf['Name'].strip().replace(' ', '').replace(',', '')
        evdat = evinf['StartTime'].strftime('%d%b%Y')
        fulldir = os.path.join(newifname, '%s_%s_%s' % (e1, evname, evdat))

        # Connect to database
        connection = None
        connection = lite.connect(newdbname)
        cursor = connection.cursor()

        # Get max gis_id currently existing
        idlist = cursor.execute("""SELECT gis_id FROM gisfiles""")
        gisids = idlist.fetchall()
        gisids = [id1[0] for id1 in gisids]
        currentgid = np.max(gisids) + 1

        # Go through and extract and rename all of the associated files, make summary file on each page
        cursor_output = cursor.execute("""SELECT gis_id FROM gisfiles WHERE event_id=?""", (e1, ))
        if cursor_output.fetchone() is not None:
            cursor = connection.cursor()
            cursor_output = cursor.execute("""SELECT gis_id, type, source, description, file_extension, apply_also FROM gisfiles WHERE event_id=?""", (e1, ))
            gisdir = os.path.join(fulldir, 'maps_gis')
            try:
                os.mkdir(gisdir)
            except Exception as e:
                print(e)
            csvfilename = os.path.join(fulldir, '%s_%s_%s_maps_gis.csv' % (e1, evname, evdat))
            if os.path.isfile(csvfilename):
                perm = 'a'
            else:
                perm = 'wb'
            with open(csvfilename, perm) as csvfile:
                num = 1
                # Put something in for the first line
                evdat = evinf['StartTime'].strftime('%d%b%Y')
                writer = csv.writer(csvfile)
                if perm == 'wb':
                    writer.writerow(['GIS files and maps of %s %s %s (Event_id %d)' % (evinf['StartTime'].strftime('%d%b%YT%H:%M'), evinf['Name'], evinf['Type'], evinf['Eid'])])
                    writer.writerow(['Note: Some of these files may be located in the folders for other events that share background information'])
                    writer.writerow(['gis_id', 'type', 'source', 'description', 'file_extension', 'apply also'])
                for temp in cursor_output:
                    # go through each entry, extract names, rename, copy over, write it down in csv
                    gis_id, type1, source, description, file_extension, apply_also = temp
                    cnt = 0
                    # get list of files if there is more than one
                    if ';' in file_extension:
                        manyfiles = file_extension.replace('; ', ';').split(';')
                        filenames = [os.path.join(relpath, mf) for mf in manyfiles]
                    else:
                        filenames = [os.path.join(relpath, file_extension)]

                    for filen in filenames:
                        fn, ext = os.path.splitext(filen)
                        indx = len(str(e1))
                        if str(e1) != os.path.basename(filen)[:indx]:
                            if 'Without.' in filen:
                                descr = 'NoOutline'
                            elif 'With' in filen:
                                descr = 'WithOutline'
                            else:
                                descr = ''
                            descr = type1.replace(' ', '').split('-')[0] + '_' + descr
                            temp = '%s_%s_%s_%s' % (e1, evname, evdat, descr)
                            if temp[-1] == '_':
                                temp = temp[:-1]
                            newfilen = os.path.join(fulldir, 'maps_gis', temp + ext.lower())
                            newfilenrel = os.path.join(shortfilen, '%s_%s_%s' % (e1, evname, evdat), 'maps_gis', temp + ext.lower())
                            num += 1
                        else:
                            temp = os.path.basename(filen)
                            newfilen = os.path.join(fulldir, 'maps_gis', temp)
                            newfilenrel = os.path.join(shortfilen, '%s_%s_%s' % (e1, evname, evdat), 'maps_gis', temp)
                        if filen != newfilen:
                            copy(filen, newfilen)
                            if cnt == 0:
                                # replace original entry
                                try:
                                    with connection:
                                        connection.execute('UPDATE gisfiles SET file_extension=? WHERE gis_id = ?', (newfilenrel, gis_id))
                                    writer.writerow([gis_id, type1, source, description, newfilenrel])
                                    cnt += 1
                                except Exception as e:
                                    print(e)
                            else:
                                # add new entry
                                try:
                                    with connection:
                                        connection.execute('INSERT INTO gisfiles(gis_id, event_id, type, source, description, file_extension, apply_also) VALUES (?,?,?,?,?,?,?)',
                                                           (currentgid, e1, type1, source, description, newfilenrel, apply_also))
                                    writer.writerow([gis_id, type1, source, description, newfilenrel])
                                except Exception as e:
                                    print(e)
                                currentgid += 1

                            if apply_also is not None:
                                aids = apply_also.replace(' ', '').split(',')
                                aids = [int(a1) for a1 in aids]
                                for a1 in aids:
                                    # Try to make file, if not already there
                                    aevinf = findsta.getEventInfo(a1, database=newdbname)
                                    aevname = aevinf['Name'].strip().replace(' ', '').replace(',', '')
                                    aevdat = aevinf['StartTime'].strftime('%d%b%Y')
                                    afulldir = os.path.join(newifname, '%s_%s_%s' % (a1, aevname, aevdat))
                                    acsvfilename = os.path.join(afulldir, '%s_%s_%s_maps_gis.csv' % (a1, aevname, aevdat))
                                    if os.path.isfile(acsvfilename):
                                        perm = 'a'
                                    else:
                                        perm = 'wb'
                                    with open(acsvfilename, perm) as acsvfile:
                                        awriter = csv.writer(acsvfile)
                                        if perm == 'wb':
                                            # if new Put something in for the first line
                                            awriter.writerow(['GIS files and maps of %s %s %s (Event_id %d)' % (aevinf['StartTime'].strftime('%d%b%YT%H:%M'), aevinf['Name'], aevinf['Type'], aevinf['Eid'])])
                                            awriter.writerow(['gis_id', 'type', 'source', 'description', 'file_extension'])
                                            awriter.writerow(['Note: Some of these files may be located in the folders for other events that share background information'])
                                        awriter.writerow([gis_id, type1, source, description, newfilenrel])

    connection.close()
    time.sleep(2.0)


def clean_information(eids, relpath, newdbname, newifname, shortfilen):
    """Make flat file of all information entries for each event and list of references, copy reference files to a top
        level references directory.

    Args:
        eids (list): list of event ids
        relpath (str): directory where infofiles is located
        newdbname (str): name of new database
        newifname (str): name of new directory for info files
        shortfilen (str): relative file path to append to new file locations when they are replaced in database

    """

    # Connect to database
    connection = None
    connection = lite.connect(newdbname)
    cursor = connection.cursor()

    # Get max ref id currently existing
    idlist = cursor.execute("""SELECT Refid FROM [reference]""")
    phoids = idlist.fetchall()
    phoids = [id1[0] for id1 in phoids]
    #currentrid = np.max(phoids) + 1

    # Go through each event
    for e1 in eids:
        # Get basic event information
        evinf = findsta.getEventInfo(e1, database=newdbname)
        evname = evinf['Name'].strip().replace(' ', '').replace(',', '')
        evdat = evinf['StartTime'].strftime('%d%b%Y')
        fulldir = os.path.join(newifname, '%s_%s_%s' % (e1, evname, evdat))

        # information table - instead of reference id, put short ref and then put those references in references csv
        # See if event has any information entries
        cursor_output = cursor.execute("""SELECT iid FROM information WHERE event_id=?""", (e1,))
        if cursor_output.fetchone() is not None:
            cursor_output = cursor.execute("""SELECT iid, ref_id, note, pho_id, GIS_id, apply_also FROM information
                                           WHERE event_id=?""", (e1,))
            data = cursor_output.fetchall()
            csvfilename = os.path.join(fulldir, '%s_%s_%s_information.csv' % (e1, evname, evdat))
            if os.path.isfile(csvfilename):
                perm = 'a'
            else:
                perm = 'wb'

            with open(csvfilename, perm) as csvfile:
                filelist = ''
                # Put something in for the first line
                writer = csv.writer(csvfile)
                if perm == 'wb':
                    writer.writerow(['Notes for %s %s %s (Event_id %d)' % (evinf['StartTime'].strftime('%d%b%YT%H:%M'),
                                    evinf['Name'], evinf['Type'], evinf['Eid'])])
                    writer.writerow(['iid', 'note', 'source', 'ref_id', 'photo_id', 'file_extension', 'url'])

                    for temp in data:
                        # go through each entry, write it down in csv
                        iid, ref_id, note, pho_id, GIS_id, apply_also = temp
                        filelist = ''
                        if GIS_id is not None:
                            # Get file location
                            gisinfo = cursor.execute("""SELECT file_extension FROM gisfiles WHERE event_id=?""", (str(e1),))
                            file_extension = gisinfo.fetchone()[0]
                            newfilenrel = os.path.join(shortfilen, 'maps_gis', '%s_%s_%s' % (e1, evname, evdat), os.path.basename(file_extension))
                            filelist = '; '.join([filelist, newfilenrel])
                        if ref_id is not None:
                            refinfo = cursor.execute("""SELECT Refid, short_ref, url FROM [reference] WHERE Refid=?""", (ref_id,))
                            ri = refinfo.fetchone()
                            rid, short_ref1, url = ri
                        else:
                            short_ref1 = ''
                            url = ''

                        if pho_id is not None:
                            temppinfo = cursor.execute("""SELECT photographer FROM photos WHERE photo_id=?
                                                       and event_id=?""", (pho_id, e1))
                            photographer = temppinfo.fetchone()[0]
                            evdat = evinf['StartTime'].strftime('%d%b%Y')
                            photoshort = photographer.replace(' ', '').replace(',', '-').strip('.')
                            temp = '%s_%s_%s_source-%s*' % (e1, evname, evdat, photoshort)
                            newfilenrel = os.path.join(shortfilen, '%s_%s_%s' % (e1, evname, evdat), 'photos_figures', temp)
                            filelist = '; '.join([filelist, newfilenrel])
                            short_ref1 = '; '.join([short_ref1, photographer])
                        if len(filelist) > 0:
                            if filelist[0] == ';':
                                filelist = filelist[2:]
                        if len(short_ref1) > 0:
                            if short_ref1[0] == ';':
                                short_ref1 = short_ref1[2:]

                        writer.writerow([iid, note.encode('utf-8'), short_ref1, ref_id, pho_id, filelist, url])

                        if apply_also is not None:
                            aids = apply_also.replace(' ', '').split(',')
                            aids = [int(a1) for a1 in aids]
                            for a1 in aids:
                                # Try to make file, if not already there
                                aevinf = findsta.getEventInfo(a1, database=newdbname)
                                aevname = aevinf['Name'].strip().replace(' ', '').replace(',', '')
                                aevdat = aevinf['StartTime'].strftime('%d%b%Y')
                                afulldir = os.path.join(newifname, '%s_%s_%s' % (a1, aevname, aevdat))
                                acsvfilename = os.path.join(afulldir, '%s_%s_%s_information.csv' % (a1, aevname, aevdat))

                                if os.path.isfile(acsvfilename):
                                    perm = 'a'
                                else:
                                    perm = 'wb'

                                with open(acsvfilename, perm) as acsvfile:
                                    awriter = csv.writer(acsvfile)
                                    if perm == 'wb':
                                        # if new Put something in for the first line
                                        awriter.writerow(['Notes for %s %s %s (Event_id %d)' % (aevinf['StartTime'].strftime('%d%b%YT%H:%M'), aevinf['Name'], aevinf['Type'], aevinf['Eid'])])
                                        awriter.writerow(['iid', 'note', 'source', 'ref_id', 'photo_id', 'file_extension', 'url'])
                                    awriter.writerow([iid, note.encode('utf-8'), short_ref1, ref_id, pho_id, filelist, url])
                        filelist = ''
                        short_ref1 = None
                        pho_id = None
                        ref_id = None
                        temp = None

    connection.close()


def clean_references(eids, newdbname, newifname, relpath):
    """Copy all references to one big file, and post to individual event files

    Args:
        eids (list): list of event ids
        newdbname (str): name of new database
        newifname (str): name of new directory for info files
        relpath (str): directory where infofiles is located

    """
    # Connect to database
    connection = None
    connection = lite.connect(newdbname)
    cursor = connection.cursor()

    cursor_output = connection.execute("""SELECT Refid FROM [reference] ORDER by short_ref""")
    ref_id1 = cursor_output.fetchall()
    rids = [e2[0] for e2 in ref_id1]
    currentrid = np.max(rids) + 1

    with open(os.path.join(os.path.dirname(newdbname), 'references.csv'), 'wb') as csvfile:
        writer2 = csv.writer(csvfile)
        # writer2.writerow(['List of all references cited'])
        writer2.writerow(['short_ref', 'Refid', 'long_ref', 'file_extension', 'url', 'event_ids'])
        for ref_id in rids:
            refinfo = cursor.execute("""SELECT event_id, apply_also, short_ref, long_ref, file_extension, url FROM [reference] WHERE Refid=?""", (ref_id,))
            ri = refinfo.fetchone()
            event_id, apply_also, short_ref1, long_ref1, file_extension1, url1 = ri
            shortsref = short_ref1.replace('(', '').replace(')', '').replace(' ', '_').replace(',', '').replace('.', '')
            if event_id is not None:
                eventlist = [int(event_id)]
                if apply_also is not None:
                    appl = apply_also.split(',')
                    for app in appl:
                        eventlist.append(int(app))
                    eventlist = sorted(eventlist)
                evstring = ', '.join([str(k) for k in eventlist])
            else:
                evstring = ''

            if file_extension1 is not None:
                cnt = 0
                # get list of files if there is more than one
                if ';' in file_extension1:
                    manyfiles = file_extension1.replace('; ', ';').split(';')
                    filenames = [os.path.join(relpath, mf) for mf in manyfiles]
                else:
                    filenames = [os.path.join(relpath, file_extension1)]
                for filen in filenames:
                    fn, ext = os.path.splitext(filen)
                    temp = '%s%s' % (shortsref, ext)
                    newfilen = os.path.join(os.path.dirname(newdbname), 'references', temp)
                    newfilenrel = os.path.join('references', temp)
                    if filen != newfilen:
                        copy(filen, newfilen)
                        if cnt == 0:
                            # replace original entry
                            #try:
                            with connection:
                                connection.execute('UPDATE [reference] SET file_extension=? WHERE Refid = ?',
                                                   (newfilenrel, ref_id))
                            writer2.writerow([short_ref1, ref_id, long_ref1.encode('utf-8'), newfilenrel, url1, evstring])
                            cnt += 1
                            #except Exception as e:
                            #    print e
                        else:
                            # add new entry
                            #try:
                            with connection:
                                connection.execute('INSERT INTO [reference](Refid, short_ref, long_ref, file_extension, url) VALUES (?,?,?,?,?)', (currentrid, short_ref1, long_ref1, newfilenrel, url1))
                            writer2.writerow([short_ref1, ref_id, long_ref1.encode('utf-8'), newfilenrel, url1, evstring])
                            cnt += 1
                            #except Exception as e:
                            #    print e
                            currentrid += 1
            else:
                writer2.writerow([short_ref1, ref_id, long_ref1.encode('utf-8'), '', url1, evstring])

    # Then do individual files
    for e1 in eids:
        # Get basic event information
        evinf = findsta.getEventInfo(e1, database=newdbname)
        evname = evinf['Name'].strip().replace(' ', '').replace(',', '')
        evdat = evinf['StartTime'].strftime('%d%b%Y')
        fulldir = os.path.join(newifname, '%s_%s_%s' % (e1, evname, evdat))

        # information table - instead of reference id, put short ref and then put those references in references csv
        # See if event has any information entries
        cursor_output = cursor.execute("""SELECT Refid FROM reference WHERE event_id=?""", (e1,))
        if cursor_output.fetchone() is not None:
            cursor_output = cursor.execute("""SELECT Refid, apply_also, short_ref, long_ref, file_extension, url FROM [reference]
                                           WHERE event_id=?""", (e1,))
            data = cursor_output.fetchall()
            csvfilename = os.path.join(fulldir, '%s_%s_%s_references.csv' % (e1, evname, evdat))
            if os.path.isfile(csvfilename):
                perm = 'a'
            else:
                perm = 'wb'

            with open(csvfilename, perm) as csvfile:
                # Put something in for the first line
                writer = csv.writer(csvfile)
                if perm == 'wb':
                    writer.writerow(['References relevant to %s %s %s (Event_id %d)' % (evinf['StartTime'].strftime('%d%b%YT%H:%M'),
                                    evinf['Name'], evinf['Type'], evinf['Eid'])])
                    writer.writerow(['short_ref', 'Refid', 'long_ref', 'file_extension', 'url'])
                    for temp in data:
                        # go through each entry, write it down in csv
                        Refid, apply_also, short_ref1, long_ref1, file_extension1, url1 = temp

                        writer.writerow([short_ref1, Refid, long_ref1.encode('utf-8'), file_extension1, url1])

                        if apply_also is not None:
                            aids = apply_also.replace(' ', '').split(',')
                            aids = [int(a1) for a1 in aids]
                            for a1 in aids:
                                aevinf = findsta.getEventInfo(a1, database=newdbname)
                                aevname = aevinf['Name'].strip().replace(' ', '').replace(',', '')
                                aevdat = aevinf['StartTime'].strftime('%d%b%Y')
                                afulldir = os.path.join(newifname, '%s_%s_%s' % (a1, aevname, aevdat))
                                acsvfilename = os.path.join(afulldir, '%s_%s_%s_references.csv' % (a1, aevname, aevdat))
                                if os.path.isfile(acsvfilename):
                                    perm = 'a'
                                else:
                                    perm = 'wb'
                                with open(acsvfilename, perm) as acsvfile:
                                    awriter = csv.writer(acsvfile)
                                    if perm == 'wb':
                                        # if new Put something in for the first line
                                        awriter.writerow(['References relevant to %s %s %s (Event_id %d)' % (aevinf['StartTime'].strftime('%d%b%YT%H:%M'), aevinf['Name'], aevinf['Type'], aevinf['Eid'])])
                                        awriter.writerow(['short_ref', 'Refid', 'long_ref', 'file_extension', 'url'])
                                    awriter.writerow([short_ref1, Refid, long_ref1.encode('utf-8'), file_extension1, url1])

    connection.close()


def clean_seismic(eids, newdbname, newifname, relpath, shortfilen):
    """Make flat file summarizing seismic detection information for each event

    Args:
        eids (list): list of event ids
        newdbname (str): name of new database
        newifname (str): name of new directory for info files
        relpath (str): directory where infofiles is located
        shortfilen (str): relative file path to append to new file locations when they are replaced in database

    """

    # Connect to database
    connection = None
    connection = lite.connect(newdbname)
    cursor = connection.cursor()

    # Go through each event
    for e1 in eids:
        # Get basic event information
        evinf = findsta.getEventInfo(e1, database=newdbname)
        evname = evinf['Name'].strip().replace(' ', '').replace(',', '')
        evdat = evinf['StartTime'].strftime('%d%b%Y')
        fulldir = os.path.join(newifname, '%s_%s_%s' % (e1, evname, evdat))

        stad = findsta.getStaInfo(e1)

        # get max distance examined -
        # searchdists = np.array([st['stasource_radius_km'] for junk, st in stad.items() if st['detect_HF'] is not None or st['detect_LP'] is not None])
        numHF = len([st for junk, st in list(stad.items()) if st['detect_HF'] == 1])
        numLP = len([st for junk, st in list(stad.items()) if st['detect_LP'] == 1])
        # don't include null, but do include nos
        # if len(searchdists) == 0:
        #     continue
        # else:
        #     maxsearch = searchdists.max()

        # seismic detection table
        cursor_output = cursor.execute("""SELECT Name, Network, Channel, LocationCode, Latitude, Longitude, Elevation_masl,
                                       stasource_radius_km, az, baz, detect_HF, detect_LP, source
                                       FROM sta_nearby, stations ON sta_nearby.station_id = stations.Sid
                                       WHERE event_id = ?
                                       ORDER BY stasource_radius_km""", (e1,))
        data = cursor_output.fetchall()
        if data is not None:
            csvfilename = os.path.join(fulldir, '%s_%s_%s_seismic_detections.csv' % (e1, evname, evdat))
            with open(csvfilename, 'wb') as csvfile:
                # Put something in for the first line
                evdat = evinf['StartTime'].strftime('%d%b%Y')
                writer = csv.writer(csvfile)
                writer.writerow(['Seismic detection for %s %s %s (Event_id %d)' % (evinf['StartTime'].strftime('%d%b%YT%H:%M'), evinf['Name'], evinf['Type'], evinf['Eid'],)])
                # writer.writerow(['Maximum distance examined: %0.1f km' % maxsearch])
                writer.writerow(['Number of High Frequency detections: %d' % numHF])
                writer.writerow(['Number of Long Period detections: %d' % numLP])
                writer.writerow(['Notes: Estimated start and end times for individual stations are only picked on stations with clear detections and valid station correction information.'])
                writer.writerow(['Detections were based on if the signal was visibly above the noise level in the frequency band specified. Detections may vary depending on the analyst and if other frequency bands are examined.'])
                writer.writerow(['Long period records were examined in displacement records, while high frequency detections were examined in velocity records.'])
                writer.writerow(['Station', 'Network', 'Channel', 'Location Code', 'Latitude', 'Longitude', 'Elevation (m a.s.l)', 'Distance (km)', 'Azimuth (deg)',
                                 'Backazimuth (deg)', 'High Frequency (HF) Detection (1-5 Hz)', 'Long Period (LP) Detection (20-60 sec)', 'data source'])
                for dat in data:
                    Name, Network, Channel, LocationCode, Latitude, Longitude, Elevation_masl, stasource_radius_km, az, baz, detect_HF, detect_LP, source = dat
                    LocationCode = str(LocationCode).zfill(2)
                    if detect_HF is not None or detect_LP is not None:
                        # if (detect_HF is None and detect_LP == 0) or (detect_HF == 0 and detect_LP is None):
                        #     continue
                        # else:
                        if detect_HF == 1:
                            detect_HF = True
                        if detect_HF == 0:
                            detect_HF = False
                        if detect_LP == 1:
                            detect_LP = True
                        if detect_LP == 0:
                            detect_LP = False
                        writer.writerow([Name, Network, Channel, LocationCode, Latitude, Longitude, Elevation_masl, stasource_radius_km, az, baz, detect_HF, detect_LP, source])
            # Copy over seismic data if necessary
            datloc = evinf['DatLocation']
            datl = datloc.split(',')
            replacewith = ''
            flag = 0
            for s1 in datl:
                if 'zip:' in s1:
                    filen = s1.split('zip:')[1].replace(' ', '')
                    fullfilename = os.path.join(relpath, filen)
                    fn, ext = os.path.splitext(fullfilename)
                    temp = '%s_%s_%s_seismic_data' % (e1, evname, evdat)
                    newfilen = os.path.join(fulldir, temp + ext.lower())
                    newfilenrel = os.path.join(shortfilen, '%s_%s_%s' % (e1, evname, evdat), temp + ext.lower())
                    copy(fullfilename, newfilen)
                    replacewith = '%s, %s' % (replacewith, newfilenrel)
                    flag = 1
                elif 'sac:' not in s1:
                    replacewith = '%s, %s' % (replacewith, s1.strip())
                    flag = 1
            if flag == 1:
                if replacewith[0] == ',':
                    replacewith = replacewith[1:]
                with connection:
                    connection.execute('UPDATE events SET DatLocation=? WHERE Eid = ?', (replacewith, e1))


def make_eventsummary(eids, newdbname, newifname, shortfilen):
    """Make single event table summary flat file

    Args:
        eids (list): list of event ids
        newdbname (str): name of new database
        newifname (str): name of new directory for info files
        shortfilen (str): relative file path to append to new file locations when they are replaced in database

    """
    # Connect to database
    connection = None
    connection = lite.connect(newdbname)
    cursor = connection.cursor()

    with open(os.path.join(os.path.dirname(newdbname), 'Events.csv'), 'wb') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['StartTime', 'EndTime', 'Event_id', 'Name', 'StateProvince', 'Country', 'Type', 'Latitude', 'Longitude', 'LocUncert_km', 'Crown_lat', 'Crown_lon', 'Tip_lat', 'Tip_lon', 'Area_total', 'Area_source', 'Area_source_low', 'Area_source_high', 'Volume', 'Volume_low', 'Volume_high', 'Mass', 'Mass_low', 'Mass_high', 'H', 'H_low', 'H_high', 'L', 'L_low', 'L_high', 'OtherDataQuality', 'LPpotential', 'maxdistHF_km', 'maxdistHF_reached', 'numHFdet', 'maxdistLP_km', 'maxdistLP_reached', 'numLPdet', 'maxdist_examined', 'DatLocation', 'FolderName'])
        writer.writerow(['UTC', 'UTC', 'integer', 'text', 'text', 'text', 'text', 'decimal degrees', 'decimal degrees', 'km', 'decimal degrees', 'decimal degrees', 'decimal degrees', 'decimal degrees', 'm^2', 'm^2', 'm^2', 'm^2', 'm^3', 'm^3', 'm^3', 'kg', 'kg', 'kg', 'm', 'm', 'm', 'm', 'm', 'm', 'integer', 'integer', 'km', 'boolean', 'integer', 'km', 'boolean', 'integer', 'km', 'text', 'text'])
        cursor_output = cursor.execute("""SELECT Eid, Name, Type, Starttime, Endtime, Latitude, Longitude, LocUncert_km, Crown_lat, Crown_lon, Tip_lat, Tip_lon, Area_total, Area_source, Area_source_low, Area_source_high, Volume, Volume_low, Volume_high, Mass, Mass_low, Mass_high, H, H_low, H_high, L, L_low, L_high, OtherDataQuality, LPpotential, maxdistHF_km, maxdistHF_reached, maxdistLP_km, maxdistLP_reached, DatLocation FROM events ORDER BY Starttime DESC""")
        data = cursor_output.fetchall()
        for dat in data:
            Eid, Name, Type, Starttime, Endtime, Latitude, Longitude, LocUncert_km, Crown_lat, Crown_lon, Tip_lat, Tip_lon, Area_total, Area_source, Area_source_low, Area_source_high, Volume, Volume_low, Volume_high, Mass, Mass_low, Mass_high, H, H_low, H_high, L, L_low, L_high, OtherDataQuality, LPpotential, maxdistHF_km, maxdistHF_reached, maxdistLP_km, maxdistLP_reached, DatLocation = dat

            # Get basic event information
            evinf = findsta.getEventInfo(Eid, database=newdbname)
            evname = evinf['Name'].strip().replace(' ', '').replace(',', '')
            evdat = evinf['StartTime'].strftime('%d%b%Y')
            fulldir = os.path.join(shortfilen, '%s_%s_%s' % (Eid, evname, evdat))
            stad = findsta.getStaInfo(Eid, database=newdbname)
            # get max distance examined -
            searchdists = np.array([st['stasource_radius_km'] for junk, st in list(stad.items()) if st['detect_HF'] is not None or st['detect_LP'] is not None])
            results = rg.search((Latitude, Longitude))
            State = results[0]['admin1']
            Country = results[0]['cc']
            if len(searchdists) > 1:
                searchdists = searchdists.max()
            else:
                searchdists = 0.

            numHF = len([st for junk, st in list(stad.items()) if st['detect_HF'] == 1])
            numLP = len([st for junk, st in list(stad.items()) if st['detect_LP'] == 1])
            writer.writerow([Starttime, Endtime, Eid, Name, State, Country, Type, Latitude, Longitude, LocUncert_km, Crown_lat, Crown_lon, Tip_lat, Tip_lon, Area_total, Area_source, Area_source_low, Area_source_high, Volume, Volume_low, Volume_high, Mass, Mass_low, Mass_high, H, H_low, H_high, L, L_low, L_high, OtherDataQuality, LPpotential, maxdistHF_reached, maxdistHF_km, numHF, maxdistLP_km, maxdistLP_reached, numLP, searchdists, DatLocation, fulldir])

    connection.close()


def main(newdbname=None, newifname=None, database='/Users/kallstadt/LSseis/landslideDatabase/lsseis.db',
         infofiles='/Users/kallstadt/LSseis/landslideDatabase/InfoFiles', shortfilen=None):
    """Run all cleaning functions

    Args:
        newdbname (str): name of new database
        newifname (str): name of new directory for info files
        database (str): full file path to sqlite3 database to clean
        infofiles (str): full file path to InfoFiles location
        shortfilen (str): relative file path to append to new file locations when they are replaced in database

    """
    eids, relpath, newdbname, newifname, shortfilen = prepare_new(newdbname=newdbname, newifname=newifname, database=database, infofiles=infofiles, shortfilen=shortfilen)
    clean_photos(eids, relpath, newdbname, newifname, shortfilen)
    clean_gis(eids, relpath, newdbname, newifname, shortfilen)
    clean_information(eids, relpath, newdbname, newifname, shortfilen)
    clean_references(eids, newdbname, newifname, relpath)
    clean_seismic(eids, newdbname, newifname, relpath, shortfilen)
    make_eventsummary(eids, newdbname, newifname, shortfilen)
