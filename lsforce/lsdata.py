from obspy.geodetics import gps2dist_azimuth
from obspy.clients.fdsn import Client
import numpy as np
import warnings

KM_PER_M = 1 / 1000  # [km/m]


class LSData:
    """Class for force inversion data (essentially a souped-up Stream).

    Attributes:
        st_orig: Original input Stream `st`
        st_proc: Stream rotated into RTZ w.r.t. `source_lat`, `source_lon`
        source_lat: See below
        source_lon: See below
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

        # Rotate into RTZ
        self.st_proc = _rotate_to_rtz(self.st_proc)

        self.st_proc.sort(keys=['distance', 'channel'])


def _rotate_to_rtz(st):
    """Rotates all components of a Stream into radial–transverse–vertical.

    This function first rotates non-standard horizontals into EN, then rotates into
    RTZ. It also checks that input components labeled as east, north, and vertical
    (e.g., BHE, BHN, BHZ, etc.) have the correct orientation.

    Args:
        st: ObsPy Stream

    Returns:
        Rotated Stream
    """

    st_rot = st.copy()  # Work on a copy of the data

    # Assemble info for inventory gather
    networks = np.unique([tr.stats.network for tr in st_rot])
    stations = np.unique([tr.stats.station for tr in st_rot])
    channels = np.unique([tr.stats.channel for tr in st_rot])

    # Grab inventory for orientation info (relying only on IRIS FDSN here!)
    client = Client('IRIS')
    inv = client.get_stations(
        network=','.join(networks),
        station=','.join(stations),
        channel=','.join(channels),
        starttime=st_rot[0].stats.starttime,
        endtime=st_rot[0].stats.endtime,
        level='channel',
    )

    # Rotate on a station-by-station basis
    for station in stations:

        st_sta = st_rot.select(station=station)  # Just select Traces for this station
        components = [tr.stats.channel[-1] for tr in st_sta]
        components.sort()  # Increasing numerical, then alphabetical order

        # ENZ case (or subset of ENZ)
        if components in (
            ['E', 'N', 'Z'],
            ['E', 'N'],
            ['E', 'Z'],
            ['N', 'Z'],
            ['E'],
            ['N'],
            ['Z'],
        ):
            # Just check for expected orientation, don't do any rotation
            for component in components:

                # Get orientation
                tr_id = st_sta.select(component=component)[0].id  # Only one Trace here!
                orient = inv.get_orientation(tr_id)
                azimuth = orient['azimuth']
                dip = orient['dip']

                # Define what it means to be "bad"
                bad_e = component == 'E' and not (azimuth == 90.0 and dip == 0.0)
                bad_n = component == 'N' and not (azimuth == 0.0 and dip == 0.0)
                bad_z = component == 'Z' and not (azimuth == 0.0 and dip == -90.0)

                # Warn!
                if bad_e or bad_n or bad_z:
                    warnings.warn(f'Bad orientation for {tr_id}\n\tazimuth = '
                                  f'{azimuth}\n\tdip = {dip}')

        # 12Z and 123 cases
        elif components in (['1', '2', 'Z'], ['1', '2', '3']):
            st_sta.rotate('->ZNE', inventory=inv)  # Rotate to ENZ (in-place change!)

        # Error out since we don't know how to handle this (yet)!
        else:
            raise ValueError(f'Unable to rotate station {station}')

        # Now the data are ostensibly in correct ENZ orientation - do rotation to RTZ
        st_sta.rotate('NE->RT')

    return st_rot
