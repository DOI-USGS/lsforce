from obspy.geodetics import gps2dist_azimuth
import numpy as np

KM_PER_M = 1/1000  # [km/m]


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
