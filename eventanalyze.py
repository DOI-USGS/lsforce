#!/usr/bin/env python

"""
Functions for analyzing events in the landslide database
"""

import sqlite3 as lite
from datetime import datetime
#import urllib2
#from obspy.iris import Client
from obspy import Stream, UTCDateTime
import numpy as np
from reviewData import reviewData
import findsta


def 


    return st


def make_measurements(event_id, straw=None, buffer_sec=None, HF=True, LP=True, HFlims=(1., 5.),
                      HFoutput='VEL', LPlims=(20., 60.), LPoutput='DISP', minradius=0.,
                      maxradius=None, database='/Users/kallstadt/LSseis/landslideDatabase/lsseis.db'):
    pass
