#!/usr/bin/env python

"""
Functions for cleaning up the database, renaming files and so on
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

# SET NULL ANY EH channels for detectLP and set NULL for any greater than maxdist

def get_exif(imgfilename):
    """Returns info from exif data of an Image, if it exists.
    Modified from https://gist.github.com/erans/983821#file-get_lat_lon_exif_pil-py-L4
    :param imgfilename: path to image filenames
    :returns:
        latitude: latitude of where photo was taken from, in degrees, from GPSTAGS
        longitude: longitude of where photo was taken from, in degrees, from GPSTAGS
        date: formatted string of date and time image was taken
    """
    try:
        image = Image.open(imgfilename)
        exif_data = {}
        info = image._getexif()
        if info:
            for tag, value in info.items():
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
    """
    Set up to make a cleaned copy of the database along with flat files for each event
    :param newdbname: path and file name for new copy of database, if None, will use current name (database) with date appended
    :param newifname: directory where new copies of info files should be placed, if None, will use current name (infofiles) with date appended
    :param database: path to current database file
    :param infofiles: directory where current info files are located
    :param shortfilen: relative file path to append to new file locations when they are replaced in database, optional
    :returns:
        eids: list of event ids
        relpath: directory where infofiles is located
        newdbname: name of new database
        newifname: name of new directory for info files
        shortfilen: relative file path to append to new file locations when they are replaced in database
    """
    relpath = os.path.dirname(infofiles)
    # Create new copy of database and new copy of InfoFiles
    if newdbname is None:
        time1 = datetime.datetime.utcnow().strftime('%d%b%YT%H%M')
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

    if newifname is None:
        time1 = datetime.datetime.utcnow().strftime('%d%b%YT%H%M')
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
        fulldir = os.path.join(newifname, '%s_%s%s' % (e1, evname, evdat))
        try:
            os.mkdir(fulldir)
        except Exception as e:
            print(e)

    return eids, relpath, newdbname, newifname, shortfilen


def get_eids(newdbname):
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
    """
    Make clean copies of all photos/figures and give them uniform names, modify database to reflect changes, create flat file for each event summarizing
    :param eids: list of event ids
    :param relpath: directory where infofiles is located
    :param newdbname: name of new database
    :param newifname: name of new directory for info files
    :param shortfilen: relative file path to append to new file locations when they are replaced in database
    """

    # Go through each event
    for e1 in eids:
        # Get basic event information
        evinf = findsta.getEventInfo(e1, database=newdbname)
        evname = evinf['Name'].strip().replace(' ', '').replace(',', '')
        evdat = evinf['StartTime'].strftime('%d%b%Y')
        fulldir = os.path.join(newifname, '%s_%s%s' % (e1, evname, evdat))

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

            csvfilename = os.path.join(fulldir, '%s_%s%s_photos_figures.csv' % (e1, evname, evdat))
            if os.path.isfile(csvfilename):
                perm = 'a'
            else:
                perm = 'wb'
            with open(csvfilename, perm) as csvfile:
                num = 1
                # Put something in for the first line
                writer = csv.writer(csvfile)
                if perm == 'wb':
                    writer.writerow(['Photos of %s %s %s (Event id %d)' % (evinf['StartTime'].strftime('%d%b%YT%H:%M'), evname, evinf['Type'], evinf['Eid'])])
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
                            temp = '%s_%s%s_source-%s_%d' % (e1, evname, evdat, photoshort, num)
                            newfilen = os.path.join(fulldir, 'photos_figures', temp + ext.lower())
                            newfilenrel = os.path.join(shortfilen, '%s_%s%s' % (e1, evname, evdat), 'photos_figures', temp + ext.lower())
                            copy(filen, newfilen)
                            if cnt == 1:
                                # replace original entry
                                try:
                                    with connection:
                                        connection.execute('UPDATE photos SET latitude=?, longitude=?, date = ?,  file_extension=? WHERE photo_id = ?', (latitude, longitude, date, newfilenrel, photo_id))
                                    writer.writerow([photo_id, photographer, description, date, latitude, longitude, newfilenrel])
                                except Exception as e:
                                    print e

                            else:
                                # add new entry
                                try:
                                    with connection:
                                        connection.execute('INSERT INTO photos(photo_id, event_id, photographer, permission, description, latitude, longitude, date, file_extension, apply_also) VALUES (?,?,?,?,?,?,?,?,?,?)',
                                                           (currentpid, e1, photographer, 1, description, latitude, longitude, date, newfilenrel, apply_also))
                                    writer.writerow([currentpid, photographer, description, date, latitude, longitude, newfilenrel])
                                except Exception as e:
                                    print e
                                currentpid += 1
                            num += 1
                            cnt += 1
                    else:
                        latitude, longitude, date = get_exif(fullpath)
                        fn, ext = os.path.splitext(fullpath)
                        shortfilen = newifname.split('landslideDatabase/')[1]
                        temp = '%s_%s%s_source-%s_%d' % (e1, evname, evdat, photoshort, num)
                        newfilen = os.path.join(fulldir, 'photos_figures', temp + ext.lower())
                        newfilenrel = os.path.join(shortfilen, '%s_%s%s' % (e1, evname, evdat), 'photos_figures', temp + ext.lower())
                        copy(fullpath, newfilen)
                        num += 1
                        # replace original entry
                        try:
                            with connection:
                                connection.execute('UPDATE photos SET latitude=?, longitude=?, date = ?,  file_extension=? WHERE photo_id = ?', (latitude, longitude, date, newfilenrel, photo_id))
                        except Exception as e:
                            print e

                    if apply_also is not None:
                        aids = apply_also.replace(' ', '').split(',')
                        aids = [int(a1) for a1 in aids]
                        for a1 in aids:
                            # Try to make file, if not already there
                            aevinf = findsta.getEventInfo(a1, database=newdbname)
                            aevname = aevinf['Name'].strip().replace(' ', '').replace(',', '')
                            aevdat = aevinf['StartTime'].strftime('%d%b%Y')
                            afulldir = os.path.join(newifname, '%s_%s%s' % (a1, aevname, aevdat))
                            acsvfilename = os.path.join(afulldir, '%s_%s%s_photos_figures.csv' % (a1, aevname, aevdat))
                            if os.path.isfile(acsvfilename):
                                perm = 'a'
                            else:
                                perm = 'wb'
                            with open(acsvfilename, perm) as acsvfile:
                                awriter = csv.writer(acsvfile)
                                if perm == 'wb':
                                    # if new Put something in for the first line
                                    awriter.writerow(['Photos of %s %s %s (Event id %d)' % (aevinf['StartTime'].strftime('%d%b%YT%H:%M'), aevname, aevinf['Type'], aevinf['Eid'])])
                                    awriter.writerow(['Note: Some of these files may be located in the folders for other events that share background information'])
                                    awriter.writerow(['photo_id', 'photographer', 'description', 'date', 'latitude', 'longitude', 'file_extension'])
                                awriter.writerow([currentpid, photographer, description, date, latitude, longitude, newfilenrel])

    connection.close()
    time.sleep(2.0)


def clean_gis(eids, relpath, newdbname, newifname, shortfilen):
    """
    Make clean copies of all gis/map figures and give them uniform names, modify database to reflect changes, create flat file for each event summarizing
    :param eids: list of event ids
    :param relpath: directory where infofiles is located
    :param newdbname: name of new database
    :param newifname: name of new directory for info files
    :param shortfilen: relative file path to append to new file locations when they are replaced in database
    """
    # Go through each event
    for e1 in eids:
        # Get basic event information
        evinf = findsta.getEventInfo(e1, database=newdbname)
        evname = evinf['Name'].strip().replace(' ', '').replace(',', '')
        evdat = evinf['StartTime'].strftime('%d%b%Y')
        fulldir = os.path.join(newifname, '%s_%s%s' % (e1, evname, evdat))

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
            csvfilename = os.path.join(fulldir, '%s_%s%s_maps_gis.csv' % (e1, evname, evdat))
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
                    writer.writerow(['GIS files and maps of %s %s %s (Event id %d)' % (evinf['StartTime'].strftime('%d%b%YT%H:%M'), evname, evinf['Type'], evinf['Eid'])])
                    writer.writerow(['Note: Some of these files may be located in the folders for other events that share background information'])
                    writer.writerow(['gis_id', 'file type', 'source', 'description', 'file extension', 'apply also'])
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
                            newfilenrel = os.path.join(shortfilen, '%s_%s%s' % (e1, evname, evdat), 'maps_gis', temp + ext.lower())
                            num += 1
                        else:
                            temp = os.path.basename(filen)
                            newfilen = os.path.join(fulldir, 'maps_gis', temp)
                            newfilenrel = os.path.join(shortfilen, '%s_%s%s' % (e1, evname, evdat), 'maps_gis', temp)
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
                                    print e
                            else:
                                # add new entry
                                try:
                                    with connection:
                                        connection.execute('INSERT INTO gisfiles(gis_id, event_id, type, source, description, file_extension, apply_also) VALUES (?,?,?,?,?,?,?)',
                                                           (currentgid, e1, type1, source, description, newfilenrel, apply_also))
                                    writer.writerow([gis_id, type1, source, description, newfilenrel])
                                except Exception as e:
                                    print e
                                currentgid += 1
                            if apply_also is not None:
                                aids = apply_also.replace(' ', '').split(',')
                                aids = [int(a1) for a1 in aids]
                                for a1 in aids:
                                    # Try to make file, if not already there
                                    aevinf = findsta.getEventInfo(a1, database=newdbname)
                                    aevname = aevinf['Name'].strip().replace(' ', '').replace(',', '')
                                    aevdat = aevinf['StartTime'].strftime('%d%b%Y')
                                    afulldir = os.path.join(newifname, '%s_%s%s' % (a1, aevname, aevdat))
                                    acsvfilename = os.path.join(afulldir, '%s_%s%s_maps_gis.csv' % (a1, aevname, aevdat))
                                    if os.path.isfile(acsvfilename):
                                        perm = 'a'
                                    else:
                                        perm = 'wb'
                                    with open(acsvfilename, perm) as acsvfile:
                                        awriter = csv.writer(acsvfile)
                                        if perm == 'wb':
                                            # if new Put something in for the first line
                                            awriter.writerow(['GIS files and maps of %s %s %s (Event id %d)' % (aevinf['StartTime'].strftime('%d%b%YT%H:%M'), aevname, aevinf['Type'], aevinf['Eid'])])
                                            awriter.writerow(['gis_id', 'file type', 'source', 'description', 'file extension'])
                                            awriter.writerow(['Note: Some of these files may be located in the folders for other events that share background information'])
                                        awriter.writerow([gis_id, type1, source, description, newfilenrel])

    connection.close()
    time.sleep(2.0)


def clean_information(eids, relpath, newdbname, newifname, shortfilen):
    """
    Make flat file of all information entries for each event and list of references, copy reference files to a top level references directory
    :param eids: list of event ids
    :param relpath: directory where infofiles is located
    :param newdbname: name of new database
    :param newifname: name of new directory for info files
    :param shortfilen: relative file path to append to new file locations when they are replaced in database
    """

    # Connect to database
    connection = None
    connection = lite.connect(newdbname)
    cursor = connection.cursor()

    # Get max ref id currently existing
    idlist = cursor.execute("""SELECT Refid FROM [references]""")
    phoids = idlist.fetchall()
    phoids = [id1[0] for id1 in phoids]
    currentrid = np.max(phoids) + 1

    # Go through each event
    for e1 in eids:
        # Get basic event information
        evinf = findsta.getEventInfo(e1, database=newdbname)
        evname = evinf['Name'].strip().replace(' ', '').replace(',', '')
        evdat = evinf['StartTime'].strftime('%d%b%Y')
        fulldir = os.path.join(newifname, '%s_%s%s' % (e1, evname, evdat))

        # information table - instead of reference id, put short ref and then put those references in references csv
        # See if event has any information entries
        cursor_output = cursor.execute("""SELECT iid FROM information WHERE event_id=?""", (e1,))
        if cursor_output.fetchone() is not None:
            cursor_output = cursor.execute("""SELECT iid, ref_id, note, pho_id, GIS_id, apply_also FROM information
                                           WHERE event_id=?""", (e1,))
            data = cursor_output.fetchall()
            csvfilename = os.path.join(fulldir, '%s_%s%s_information.csv' % (e1, evname, evdat))
            if os.path.isfile(csvfilename):
                perm = 'a'
            else:
                perm = 'wb'

            refcsvfilename = os.path.join(fulldir, '%s_%s%s_references.csv' % (e1, evname, evdat))
            if os.path.isfile(refcsvfilename):
                perm2 = 'a'
            else:
                perm2 = 'wb'

            with open(csvfilename, perm) as csvfile:
                filelist = ''
                # Put something in for the first line
                writer = csv.writer(csvfile)
                if perm == 'wb':
                    writer.writerow(['Notes for %s %s %s (Event id %d)' % (evinf['StartTime'].strftime('%d%b%YT%H:%M'),
                                    evname, evinf['Type'], evinf['Eid'])])
                    writer.writerow(['iid', 'note', 'source', 'file location'])

                with open(refcsvfilename, perm2) as refcsvfile:
                    writer2 = csv.writer(refcsvfile)
                    if perm2 == 'wb':
                        writer2.writerow(['References for %s %s %s (Event id %d)' % (evinf['StartTime'].strftime('%d%b%YT%H:%M'),
                                         evname, evinf['Type'], evinf['Eid'])])
                        writer2.writerow(['Refid', 'short reference', 'full reference', 'file location', 'url'])
                    for temp in data:
                        # go through each entry, write it down in csv
                        iid, ref_id, note, pho_id, GIS_id, apply_also = temp
                        if GIS_id is not None:
                            # Get file location
                            gisinfo = cursor.execute("""SELECT file_extension FROM gisfiles WHERE event_id=?""", (str(e1),))
                            file_extension = gisinfo.fetchone()
                            newfilenrel = os.path.join(shortfilen, '%s_%s%s' % (e1, evname, evdat), os.path.basename(file_extension[0]))
                            filelist = filelist + ';' + newfilenrel
                        if ref_id is not None:
                            refinfo = cursor.execute("""SELECT Refid, short_ref, long_ref, file_extension, url FROM [references] WHERE Refid=?""", (ref_id,))
                            ri = refinfo.fetchone()
                            rid, short_ref1, long_ref1, file_extension1, url1 = ri
                            shortsref = short_ref1.replace('(', '').replace(')', '').replace(' ', '_').replace(',', '').replace('.', '')
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
                                                connection.execute('UPDATE [references] SET file_extension=? WHERE Refid = ?',
                                                                   (newfilenrel, rid))
                                            writer2.writerow([rid, short_ref1, long_ref1.encode('utf-8'), newfilenrel, url1])
                                            cnt += 1
                                            #except Exception as e:
                                            #    print e
                                        else:
                                            # add new entry
                                            #try:
                                            with connection:
                                                connection.execute('INSERT INTO [references](Refid, short_ref, Permission, long_ref, file_extension, url) VALUES (?,?,?,?,?,?)', (currentrid, short_ref1, 1, long_ref1, newfilenrel, url1))
                                            writer2.writerow([rid, short_ref1, long_ref1.encode('utf-8'), newfilenrel, url1])
                                            cnt += 1
                                            #except Exception as e:
                                            #    print e
                                            currentrid += 1
                            else:
                                writer2.writerow([rid, short_ref1, long_ref1.encode('utf-8'), '', url1])
                        else:
                            short_ref1 = ''
                        if pho_id is not None:
                            temppinfo = cursor.execute("""SELECT photographer FROM photos WHERE photo_id=?
                                                       and event_id=?""", (pho_id, e1))
                            photographer = temppinfo.fetchone()
                            evdat = evinf['StartTime'].strftime('%d%b%Y')
                            photoshort = photographer[0].replace(' ', '').replace(',', '-').strip('.')
                            temp = '%s_%s%s_source-%s*' % (e1, evname, evdat, photoshort)
                            newfilenrel = os.path.join(shortfilen, '%s_%s%s' % (e1, evname, evdat), 'photos_figures', temp)
                            filelist = filelist + ';' + newfilenrel
                        pho_id = None
                        ref_id = None
                        temp = None
                        if len(filelist) > 0:
                            if filelist[0] == ';':
                                filelist = filelist[1:]
                        writer.writerow([iid, note.encode('utf-8'), short_ref1, filelist])
                if apply_also is not None:
                    aids = apply_also.replace(' ', '').split(',')
                    aids = [int(a1) for a1 in aids]
                    for a1 in aids:
                        # Try to make file, if not already there
                        aevinf = findsta.getEventInfo(a1, database=newdbname)
                        aevname = aevinf['Name'].strip().replace(' ', '').replace(',', '')
                        aevdat = aevinf['StartTime'].strftime('%d%b%Y')
                        afulldir = os.path.join(newifname, '%s_%s%s' % (a1, aevname, aevdat))
                        acsvfilename = os.path.join(afulldir, '%s_%s%s_information.csv' % (a1, aevname, aevdat))
                        arefcsvfilename = os.path.join(afulldir, '%s_%s%s_references.csv' % (a1, aevname, aevdat))

                        if os.path.isfile(acsvfilename):
                            perm = 'a'
                        else:
                            perm = 'wb'
                        if os.path.isfile(arefcsvfilename):
                            perm2 = 'a'
                        else:
                            perm2 = 'wb'
                        with open(acsvfilename, perm) as acsvfile:
                            awriter = csv.writer(acsvfile)
                            if perm == 'wb':
                                # if new Put something in for the first line
                                awriter.writerow(['Notes for %s %s %s (Event id %d)' % (aevinf['StartTime'].strftime('%d%b%YT%H:%M'), aevname, aevinf['Type'], aevinf['Eid'])])
                                awriter.writerow(['iid', 'note', 'source', 'file location'])
                            awriter.writerow([iid, note.encode('utf-8'), short_ref1, filelist])
                            with open(arefcsvfilename, perm2) as racsvfile:
                                rwriter = csv.writer(racsvfile)
                                if perm2 == 'wb':
                                    rwriter.writerow(['References for %s %s %s (Event id %d)' % (aevinf['StartTime'].strftime('%d%b%YT%H:%M'),
                                                     aevname, aevinf['Type'], aevinf['Eid'])])
                                    rwriter.writerow(['Refid', 'short reference', 'full reference', 'file location', 'url'])
                                rwriter.writerow([rid, short_ref1, long_ref1.encode('utf-8'), '', url1])

    connection.close()


def clean_seismic(eids, newdbname, newifname, relpath, shortfilen):
    """
    Make flat file summarizing seismic detection information for each event
    :param eids: list of event ids
    :param newdbname: name of new database
    :param newifname: name of new directory for info files
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
        fulldir = os.path.join(newifname, '%s_%s%s' % (e1, evname, evdat))

        stad = findsta.getStaInfo(e1)

        # get max distance examined -
        searchdists = np.array([st['stasource_radius_km'] for junk, st in stad.items() if st['detect_HF'] is not None or st['detect_LP'] is not None])
        numHF = len([st for junk, st in stad.items() if st['detect_HF'] == 1])
        numLP = len([st for junk, st in stad.items() if st['detect_LP'] == 1])
        # don't include null, but do include nos
        if len(searchdists) == 0:
            continue
        else:
            maxsearch = searchdists.max()

        # seismic detection table
        cursor_output = cursor.execute("""SELECT Name, Network, Channel, LocationCode, Latitude, Longitude, Elevation_masl,
                                       stasource_radius_km, az, baz, detect_HF, starttimeHF, endtimeHF, absmaxampHF, detect_LP, starttimeLP, endtimeLP, absmaxampLP, source
                                       FROM sta_nearby, stations ON sta_nearby.station_id = stations.Sid
                                       WHERE event_id = ?
                                       ORDER BY stasource_radius_km""", (e1,))
        data = cursor_output.fetchall()
        if data is not None:
            csvfilename = os.path.join(fulldir, '%s_%s%s_seismic_detections.csv' % (e1, evname, evdat))
            with open(csvfilename, 'wb') as csvfile:
                # Put something in for the first line
                evdat = evinf['StartTime'].strftime('%d%b%Y')
                writer = csv.writer(csvfile)
                writer.writerow(['Seismic detection for %s %s %s (Event id %d)' % (evinf['StartTime'].strftime('%d%b%YT%H:%M'), evname, evinf['Type'], evinf['Eid'],)])
                writer.writerow(['Maximum distance examined: %0.1f km' % maxsearch])
                writer.writerow(['Number of High Frequency detections: %d' % numHF])
                writer.writerow(['Number of Long Period detections: %d' % numLP])
                writer.writerow(['Notes: Estimated start and end times for individual stations are only picked on stations with clear detections and valid station correction information.'])
                writer.writerow(['Detections and picks were made visually based on if and when the signal was above the noise level in the frequency band specified. These selections may vary depending on the analyst and may differ if other frequency bands are examined.'])
                writer.writerow(['Long period records were examined in displacement records, while high frequency detections were examined in velocity records.'])
                writer.writerow(['Station', 'Network', 'Channel', 'Location Code', 'Latitude', 'Longitude', 'Elevation (m a.s.l)', 'Distance (km)', 'Azimuth (deg)',
                                 'Backazimuth (deg)', 'High Frequency (HF) Detection (1-5 Hz)', 'Approximate HF Start time', 'Approximate HF End time', 'HF max amplitude (m/s)',
                                 'Long Period (LP) Detection (20-60 sec)', 'Approximate HF Start time', 'Approximate HF End time', 'LP max amplitude (m)', 'data source'])
                for dat in data:
                    Name, Network, Channel, LocationCode, Latitude, Longitude, Elevation_masl, stasource_radius_km, az, baz, detect_HF, starttimeHF, endtimeHF, absmaxampHF, detect_LP, starttimeLP, endtimeLP, absmaxampLP, source = dat
                    if detect_HF is not None or detect_LP is not None:
                        if (detect_HF is None and detect_LP == 0) or (detect_HF == 0 and detect_LP is None):
                            continue
                        else:
                            if detect_HF == 1:
                                detect_HF = True
                            if detect_HF == 0:
                                detect_HF = False
                            if detect_LP == 1:
                                detect_LP = True
                            if detect_LP == 0:
                                detect_LP = False
                            writer.writerow([Name, Network, Channel, LocationCode, Latitude, Longitude, Elevation_masl, stasource_radius_km, az, baz, detect_HF, starttimeHF, endtimeHF, absmaxampHF, detect_LP, starttimeLP, endtimeLP, absmaxampLP, source])
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
                    temp = '%s_%s%s_seismic_data' % (e1, evname, evdat)
                    newfilen = os.path.join(fulldir, temp + ext.lower())
                    newfilenrel = os.path.join(shortfilen, '%s_%s%s' % (e1, evname, evdat), temp + ext.lower())
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


def make_eventsummary(eids, newdbname, newifname):
    """
    Make single event table summary flat file
    :param eids: list of event ids
    :param newdbname: name of new database
    :param newifname: name of new directory for info files
    """
    # Connect to database
    connection = None
    connection = lite.connect(newdbname)
    cursor = connection.cursor()

    with open(os.path.join(os.path.dirname(newdbname), 'Events.csv'), 'wb') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Approx start time (UTC)', 'Approx end time (UTC)', 'Event id', 'Name', 'State/Province', 'Country', 'Type', 'Latitude', 'Longitude', 'Location uncertainty (km)', 'Crown latitude', 'Crown longitude', 'Tip latitude', 'Toe latitude', 'Total area (m2)', 'Source area (m2)', 'Source area lower bound (m2)', 'Source area upper bound (m2)', 'Volume (m3)', 'Volume lower bound (m3)', 'Volume upper bound (m3)', 'Mass (kg)', 'Mass lower bound (kg)', 'Mass upper bound (kg)', 'H (m)', 'H lower bound (m)', 'H upper bound (m)', 'L (m)', 'L lower bound (m)', 'L upper bound (m)', 'Available data quality', 'Long period detections?', 'High frequency (1-5 Hz) maximum observed distance (km)', 'High frequency limit reached?', 'Number of high frequency detections', 'Long period (20-60 sec) maximum observed distance (km)', 'Long period limit reached?', 'Number of long period detections', 'Maximum distance examined (km)', 'Seismic data location(s)', 'Folder name'])
        cursor_output = cursor.execute("""SELECT Eid, Name, Type, Starttime, Endtime, Latitude, Longitude, 'LocUncert_km', Crown_lat, Crown_lon, Tip_lat, Tip_lon, Area_total, Area_source, Area_source_low, Area_source_high, Volume, Volume_low, Volume_high, Mass, Mass_low, Mass_high, H, H_low, H_high, L, L_low, L_high, OtherDataQuality, LPpotential, maxdistHF_km, maxdistHF_reached, maxdistLP_km, maxdistLP_reached, DatLocation FROM events ORDER BY Starttime""")
        data = cursor_output.fetchall()
        for dat in data:
            Eid, Name, Type, Starttime, Endtime, Latitude, Longitude, LocUncert_km, Crown_lat, Crown_lon, Tip_lat, Tip_lon, Area_total, Area_source, Area_source_low, Area_source_high, Volume, Volume_low, Volume_high, Mass, Mass_low, Mass_high, H, H_low, H_high, L, L_low, L_high, OtherDataQuality, LPpotential, maxdistHF_km, maxdistHF_reached, maxdistLP_km, maxdistLP_reached, DatLocation = dat

            # Get basic event information
            evinf = findsta.getEventInfo(Eid, database=newdbname)
            evname = evinf['Name'].strip().replace(' ', '').replace(',', '')
            evdat = evinf['StartTime'].strftime('%d%b%Y')
            fulldir = os.path.join(newifname, '%s_%s%s' % (Eid, evname, evdat))
            stad = findsta.getStaInfo(Eid, database=newdbname)
            # get max distance examined -
            searchdists = np.array([st['stasource_radius_km'] for junk, st in stad.items() if st['detect_HF'] is not None or st['detect_LP'] is not None])
            results = rg.search((Latitude, Longitude))
            State = results[0]['admin1']
            Country = results[0]['cc']
            if len(searchdists) > 1:
                searchdists = searchdists.max()
            else:
                searchdists = 0.
            numHF = len([st for junk, st in stad.items() if st['detect_HF'] == 1])
            numLP = len([st for junk, st in stad.items() if st['detect_LP'] == 1])
            writer.writerow([Starttime, Endtime, Eid, Name, State, Country, Type, Latitude, Longitude, LocUncert_km, Crown_lat, Crown_lon, Tip_lat, Tip_lon, Area_total, Area_source, Area_source_low, Area_source_high, Volume, Volume_low, Volume_high, Mass, Mass_low, Mass_high, H, H_low, H_high, L, L_low, L_high, OtherDataQuality, LPpotential, maxdistHF_reached, maxdistHF_km, numHF, maxdistLP_km, maxdistLP_reached, numLP, searchdists, DatLocation, fulldir])

    connection.close()


def main(newdbname=None, newifname=None, database='/Users/kallstadt/LSseis/landslideDatabase/lsseis.db',
         infofiles='/Users/kallstadt/LSseis/landslideDatabase/InfoFiles', shortfilen=None):
    """
    Run all cleaning functions
    :param newdbname: name of new database
    :param newifname: name of new directory for info files
    :param newdbname: name of new database
    :param newifname: name of new directory for info files
    :param shortfilen: relative file path to append to new file locations when they are replaced in database
    """
    eids, relpath, newdbname, newifname, shortfilen = prepare_new(newdbname=newdbname, newifname=newifname, database=database, infofiles=infofiles, shortfilen=None)
    #clean_photos(eids, relpath, newdbname, newifname, shortfilen)
    #clean_gis(eids, relpath, newdbname, newifname, shortfilen)
    #clean_information(eids, relpath, newdbname, newifname, shortfilen)
    clean_seismic(eids, newdbname, newifname, relpath, shortfilen)
    #make_eventsummary(eids, newdbname, newifname)
