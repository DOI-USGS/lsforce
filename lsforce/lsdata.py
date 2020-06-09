from obspy.geodetics import gps2dist_azimuth
from obspy.signal.rotate import rotate2zne
from obspy import Stream, Trace
import urllib
import numpy as np

KM_PER_M = 1 / 1000  # [km/m]


class LSData:
    """Class for force inversion data (essentially a souped-up Stream).

    Attributes:
        st_orig:
        st_proc:
        source_lat:
        source_lon:
    """

    def __init__(self, st, source_lat, source_lon):
        """
        Args:
            st: ObsPy Stream object with ``tr.stats.latitude`` and
                ``tr.stats.longitude`` defined
            source_lat (int or float): Latitude in decimal degrees of centroid of
                landslide source location
            source_lon (int or float): Longitude in decimal degrees of centroid of
                landslide source location
        """

        # Before we do anything, verify that tr.stats.latitude/longitude are defined
        for tr in st:
            assert hasattr(tr.stats, 'latitude') and hasattr(tr.stats, 'longitude'),\
                f'tr.stats.latitude/longitude not present for {tr.id}'

        self.st_orig = st.copy()  # Save a copy of the original input Stream
        self.st_proc = st.copy()  # Work on a copy of the input Stream

        self.source_lat = source_lat
        self.source_lon = source_lon

        # Define distance (in km), azimuth, and back_azimuth
        for tr in self.st_proc:
            dist, az, baz = gps2dist_azimuth(
                self.source_lat, self.source_lon, tr.stats.latitude, tr.stats.longitude
            )
            tr.stats.distance = dist * KM_PER_M  # [km]
            tr.stats.azimuth = az  # [deg]
            tr.stats.back_azimuth = baz  # [deg]

        # Rotate into Z–R–T (on a station-by-station basis!)
        for station in np.unique([tr.stats.station for tr in self.st_proc]):
            self.st_proc.select(station=station).rotate('NE->RT')

        self.st_proc.sort(keys=['distance', 'channel'])


def rotate(st, back_azimuth=None):
    """
    rotate all components of st that can be rotated, to radial and transverse

    Args:
        st: obspy stream object to rotate
        back_azimuth (list): Not required if backaz already attached to st stats, list
            of backazimuths corresponding to st
    """

    # implant back_azimuth in st's
    if back_azimuth:
        for i, trace in enumerate(st):
            trace.stats.back_azimuth = back_azimuth[i]
            st[i] = trace
    else:
        try:
            st[0].stats.back_azimuth
        except:
            print('need to attach back_azimuth')
            return
    # get list of station location code pairs present
    staloc = list(
        set([trace.stats.station + '.' + trace.stats.location for trace in st])
    )
    st_rotated = Stream(traces=Trace())  # initialize, pop this one off later
    for station in staloc:
        # get all components with that station name
        loc = station.split('.')[1]
        if loc == '':
            loc = None
        st_temp = st.select(
            station=station.split('.')[0], location=loc
        ).copy()  # [trace for trace in st if station in trace.stat.station]
        if len(st_temp) == 3:  # if len 3, put z down, rotate horizontals
            try:
                z = st_temp.select(component='Z').copy()
                st_rotated = st_rotated + z.copy()
                chans = [trace.stats.channel for trace in st_temp]
                # if not BHN BHE
                if 'H1' in str(chans) and 'H2' in str(chans):
                    try:
                        for k, trace in enumerate(st_temp):
                            if trace.stats.location == '':
                                loc = '--'
                            else:
                                loc = trace.stats.location
                            url = (
                                'http://service.iris.edu/fdsnws/station/1/query?net=%s&sta=%s&loc=%s&cha=%s&level=channel&format=text&includecomments=true&nodata=404'
                                % (
                                    trace.stats.network,
                                    trace.stats.station,
                                    loc,
                                    trace.stats.channel,
                                )
                            )
                            temp = urllib.request.urlopen(url)
                            file1 = temp.read()
                            lines = [line.split('|') for line in file1.split('\n')[1:]]
                            trace.stats.cmpaz = float(lines[0][8])
                            st_temp[k] = trace

                        z1, n1, e1 = rotate2zne(
                            z[0].data,
                            0,
                            -90,
                            st_temp.select(component='1')[0].data,
                            st_temp.select(component='1')[0].stats.cmpaz,
                            0,
                            st_temp.select(component='2')[0].data,
                            st_temp.select(component='2')[0].stats.cmpaz,
                            0,
                        )
                        st_temp.select(component='1')[0].data = n1
                        st_temp.select(component='1')[0].stats.channel = 'BHN'
                        st_temp.select(component='2')[0].data = e1
                        st_temp.select(component='2')[0].stats.channel = 'BHE'
                    except:
                        print(
                            'couldnt get cmpaz orientation from IRIS, rotation failed'
                        )
                        continue
                st_h = (
                    st_temp.select(component='N').copy()
                    + st_temp.select(component='E').copy()
                )
                st_h.rotate('NE->RT')
                st_rotated = st_rotated + st_h.copy()
            except:
                print('error in rotating for ' + station + ' -skipping')
        elif len(st_temp) == 1:  # if len 1, put z down, continue
            z = st_temp.select(component='Z')
            st_rotated = st_rotated + z.copy()
        elif len(st_temp) == 2:  # if len 2, probably horizontal components
            try:
                st_h = (
                    st_temp.select(component='N').copy()
                    + st_temp.select(component='E').copy()
                )
                st_h.rotate('NE->RT')
                st_rotated = st_rotated + st_h.copy()
            except:
                print(('weird number of components for ' + station + ' -skipping'))
        else:
            print(('weird number of components for ' + station + ' -skipping'))
    st_rotated.pop(0)  # pop off the placeholder

    return st_rotated
