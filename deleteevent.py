#!/usr/bin/env python

"""
Functions for removing an event and associated entries from the database
"""

import sqlite3 as lite
import os
import shutil
import datetime


def remove_events(event_ids, savecopy=None, gisfiles=False, photos=False, information=False, database='/Users/kallstadt/LSseis/landslideDatabase/lsseis.db'):
    """
    Removes event_ids from event table and also removes any entries in sta_nearby table for this event. It also can remove entries in gisfiles, photos, or information tables if desired.
    savecopy = full file name for a copy to save of database in case something goes wrong with .db extension
    """
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
            print e
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
