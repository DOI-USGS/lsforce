import numpy as np
from cartopy import crs as ccrs
from cartopy import feature as cfeature
from matplotlib import dates as mdates
from matplotlib import patheffects as path_effects
from matplotlib import pyplot as plt
from obspy import Stream, Trace
from obspy.clients.fdsn import RoutingClient
from obspy.geodetics import gps2dist_azimuth

DETREND_POLY_ORDER = 2

KM_PER_M = 1 / 1000  # [km/m]

# Map channel sets to colors for station map
CHANNEL_COLORS = dict(
    Z='red',
    R='salmon',
    T='darkred',
    ZR='blue',
    ZT='skyblue',
    RT='darkblue',
    ZRT='green',
)


class LSData:
    r"""Class for force inversion data that is an extension of an ObsPy Stream.

    Attributes:
        st_orig (:class:`~obspy.core.stream.Stream`): Original input Stream `st`.
        st_proc (:class:`~obspy.core.stream.Stream`): Stream rotated into RTZ (radial,
            transverse, vertical) relative to `source_lat`, `source_lon`.
        source_lat (float): Latitude in decimal degrees of centroid of landslide
            source location.
        source_lon (float): Longitude in decimal degrees of centroid of landslide
            source location.
    """

    def __init__(
        self,
        st,
        source_lat,
        source_lon,
        remove_response=True,
        remove_response_kwargs=None,
        skip_zne_rotation=False,
    ):
        r"""Create an LSData object.

        Args:
            st (:class:`~obspy.core.stream.Stream`): Stream object with
                ``tr.stats.latitude`` and ``tr.stats.longitude`` defined and station
                response info attached to each trace in the Stream.
            source_lat (float): Latitude in decimal degrees of centroid of landslide
                source location
            source_lon (float): Longitude in decimal degrees of centroid of landslide
                source location
            remove_response (bool): Correct for station response to displacement units.
                Set to `False` to handle response removal manually at an earlier step.
            remove_response_kwargs (dict): Dictionary of keyword arguments to pass
                to :meth:`~obspy.core.trace.Trace.remove_response`. No effect if
                `remove_response` is `False`. Note that ObsPy's default value for the
                `water_level` kwarg is 60 dB, but this may inadvertently remove some
                long period energy that the user wants to keep, so set this parameter
                carefully. Try a higher number (e.g., 100 dB) to preserve more longer
                period energy. Setting to `None` is recommended if pre-filtering is
                applied or if fine control on the frequency limits is desired
            skip_zne_rotation (bool): If `True`, then the ->ZNE rotation step is
                skipped. This is a necessary flag if the stations used do not have
                metadata in the irisws-fedcatalog (e.g., for synthetic cases)
        """

        # Type check for source coordinates
        if type(source_lat) is not float or type(source_lon) is not float:
            raise TypeError('source_lat and source_lon must be floats.')

        # Verify that tr.stats.latitude/longitude are defined
        for tr in st:
            assert hasattr(tr.stats, 'latitude') and hasattr(tr.stats, 'longitude')

        # Ensure that if any horizontal is present, *both* are present
        for station in np.unique([tr.stats.station for tr in st]):
            st_sta = st.select(station=station)
            id_string = '\n'.join(['\t' + tr.id for tr in st_sta])
            id_msg = f'Traces involved:\n{id_string}'
            components = np.unique([tr.stats.component for tr in st_sta])
            if components.size == 1:
                msg = 'If a station has only one component, that component must be vertical!'
                assert (
                    components[0] == 'Z'
                ), f'{msg}\nComponents found: {components}\n{id_msg}'
            elif components.size == 2:
                msg = 'If a station has two components, they must both be horizontal!'
                assert (
                    'Z' not in components
                ), f'{msg}\nComponents found: {components}\n{id_msg}'

        if remove_response_kwargs is None:
            remove_response_kwargs = dict()
        remove_response_kwargs.pop('output', None)  # Remove "output" kwarg, if provided

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
        self.st_proc = _rotate_to_rtz(self.st_proc, skip_zne_rotation)

        # Now that we're rotated, sort
        self.st_proc.sort(keys=['distance', 'station', 'channel'])

        if remove_response:
            self.st_proc.detrend('polynomial', order=DETREND_POLY_ORDER)
            self.st_proc.remove_response(output='DISP', **remove_response_kwargs)

        # Set the private _is_pre_filt attribute, which is used in lsforce.py
        if remove_response and remove_response_kwargs.get('pre_filt') is not None:
            self._is_pre_filt = True
        else:
            self._is_pre_filt = False

    def plot_data(
        self, equal_scale=True, period_range=None, filter_order=2, zerophase=True
    ):
        r"""Create a record section plot of waveforms in `st_proc`.

        Args:
            equal_scale (bool): If `True`, all plots will share the same y-axis scale
            period_range (list or tuple): If not `None`, filter the data between
                `period_range[0]` and `period_range[1]`, given in seconds
            filter_order (int): Order of filter applied over period_range (no effect
                if `period_range` is `None`)
            zerophase (bool): If `True`, zero-phase filtering will be used (no effect
                if `period_range` is `None`)

        Returns:
            :class:`~matplotlib.figure.Figure`: Output figure handle
        """

        st_plot = self.st_proc.copy()  # Make a copy to manipulate

        if period_range:
            assert np.atleast_1d(period_range).size == 2, 'len(period_range) must be 2'
            period_range = sorted(period_range)  # Ensure minimum period is first
            st_plot.filter(
                'bandpass',
                freqmin=1 / period_range[1],
                freqmax=1 / period_range[0],
                corners=filter_order,
                zerophase=zerophase,
            )
            filter_string = '{}–{} s bandpass'.format(*period_range)
        else:
            filter_string = 'Unfiltered'

        if not equal_scale:
            st_plot.normalize()
            amp_string = 'normalized'
        else:
            amp_string = 'equal scale'

        # Spacing between traces is twice the max amp encountered
        spacing = np.max([np.abs(tr.data).max() for tr in st_plot]) * 2

        fig, ax = plt.subplots(figsize=(8, 12))

        offset = 0
        yticks = []
        yticklabels = []

        for tr in st_plot:
            ax.plot(
                tr.times('matplotlib'), tr.data + offset, color='black', linewidth=1
            )
            label = (
                f'{tr.stats.network}.{tr.stats.station} ({tr.stats.channel[-1]}) '
                f'– {tr.stats.distance:.1f} km'
            )
            yticklabels.append(label)
            yticks.append(offset)
            offset -= spacing

        # Misc. tweaks
        ax.set_xlim(
            st_plot[0].stats.starttime.matplotlib_date,
            st_plot[0].stats.endtime.matplotlib_date,
        )
        ax.set_ylim(yticks[-1] - spacing / 2, yticks[0] + spacing / 2)
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticklabels)
        ax.tick_params(left=False, top=True, which='both')
        ax.yaxis.set_tick_params(length=0, pad=4)

        # Time axis stuff
        locator = mdates.AutoDateLocator()
        formatter = mdates.ConciseDateFormatter(locator)
        ax.xaxis.set_major_locator(locator)
        ax.xaxis.set_major_formatter(formatter)

        ax.set_title(f'{filter_string} ({amp_string})', pad=10)

        fig.tight_layout()
        fig.show()

        return fig

    def plot_stations(self, region=None, label_stations=False, gshhs_scale='auto'):
        r"""Create a map showing stations and event location.

        Args:
            region (list or tuple): Array of the form [lonmin, lonmax, latmin, latmax]
                specifying the desired map region in decimal degrees. If `None`, we
                automatically pick a region that includes the event and stations
            label_stations (bool): If `True`, label stations with their codes
            gshhs_scale (str): Resolution for coastlines; one of `'auto'`, `'coarse'`,
                `'low'`, `'intermediate'`, `'high'`, or `'full'`

        Returns:
            :class:`~matplotlib.figure.Figure`: Output figure handle
        """

        # Automatically determine a nice region, if one not explicitly provided
        if not region:
            lons = [tr.stats.longitude for tr in self.st_proc] + [self.source_lon]
            lats = [tr.stats.latitude for tr in self.st_proc] + [self.source_lat]
            region = [
                np.floor(np.min(lons)),
                np.ceil(np.max(lons)),
                np.floor(np.min(lats)),
                np.ceil(np.max(lats)),
            ]

        proj = ccrs.AlbersEqualArea(
            central_longitude=np.mean(region[:2]),
            central_latitude=np.mean(region[2:]),
            standard_parallels=[np.min(region[2:]), np.max(region[2:])],
        )

        fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(projection=proj))

        # Add geographic context
        ax.add_feature(
            cfeature.GSHHSFeature(scale=gshhs_scale),
            facecolor=cfeature.COLORS['land'],
            zorder=1,
        )
        ax.patch.set_facecolor(cfeature.COLORS['water'])
        ax.add_feature(
            cfeature.LAKES,
            facecolor=cfeature.COLORS['water'],
            edgecolor='black',
            zorder=2,
        )

        # Universal scatter properties
        scatter_kwargs = dict(s=100, zorder=4, transform=ccrs.PlateCarree())

        # Plot event location
        ax.scatter(
            self.source_lon,
            self.source_lat,
            color='black',
            label='Event',
            **scatter_kwargs,
        )

        # Plot stations, colored by channels
        for current_channel_set, color in CHANNEL_COLORS.items():
            station_lons = []
            station_lats = []
            station_labels = []
            for station in np.unique([tr.stats.station for tr in self.st_proc]):
                station_st = self.st_proc.select(station=station)
                # If this station's set of unique channels match the channel set we're on
                station_channel_set = np.unique(
                    [tr.stats.channel[-1] for tr in station_st]
                )
                if sorted(station_channel_set) == sorted(current_channel_set):
                    station_lons.append(station_st[0].stats.longitude)
                    station_lats.append(station_st[0].stats.latitude)
                    station_labels.append(station_st[0].stats.station)
            if len(station_lons) > 0:
                # Plot and label according to CHANNEL_COLORS
                ax.scatter(
                    station_lons,
                    station_lats,
                    color=color,
                    marker='v',
                    edgecolors='black',
                    label=current_channel_set,
                    **scatter_kwargs,
                )
                if label_stations:
                    for lon, lat, label in zip(
                        station_lons, station_lats, station_labels
                    ):
                        t = ax.text(
                            lon,
                            lat,
                            '   ' + label,
                            va='center',
                            color='white',
                            transform=ccrs.PlateCarree(),
                        )
                        # Outline text
                        t.set_path_effects(
                            [
                                path_effects.Stroke(linewidth=2, foreground='black'),
                                path_effects.Normal(),
                            ]
                        )

        ax.set_extent(region, crs=ccrs.PlateCarree())
        ax.gridlines(draw_labels=True, zorder=3)
        ax.legend()

        fig.canvas.draw()
        fig.tight_layout(pad=2)
        fig.show()

        return fig


def make_lsdata_syn(
    inv,
    fake_station_dict,
    source_lat,
    source_lon,
    data_length_seconds,
    skip_zne_rotation=True,
    fake_start_time=None,
    dummy_sampling_rate=100.0,
):
    r"""Wrapper which creates an :class:`~lsforce.lsdata.LSData` object for
    forward modeling applications.

    Args:
        inv (:class:`~obspy.core.inventory.inventory.Inventory`): ObsPy Inventory
            object containing the **real** stations for which synthetics should be
            computed
        fake_station_dict (dict): Dictionary with keys specifying station names, and
            values as two-element lists [latitude, longitude], of **fake** stations
            for which synthetics should be computed
        source_lat (float): Latitude in decimal degrees of centroid of landslide
            source location
        source_lon (float): Longitude in decimal degrees of centroid of landslide
            source location
        data_length_seconds (int or float): [s] Length of synthetic data
        skip_zne_rotation (bool): Should be True if stations are entirely fake,
            but set to False if using a real station inventory and some stations
            might have BH1/BH2 channels instead of E and N.
        fake_start_time (:class:`~obspy.core.utcdatetime.UTCDateTime` or None): Fake start time to use if using
            real stations and skip_zne_rotation is set to False. This is
            needed so the right 'real' inventory data is fetched.
        dummy_sampling_rate (float): Samples per second for dummy data. If
            not defined, will use 100. Default is set artificially high, this
            determines the SEED band code of "H"
    """

    # Ensure there's something to use
    if inv is None and fake_station_dict is None:
        raise ValueError('No station information provided!')

    # Create dummy Stream object
    st = Stream()

    # Add Traces from inv, if it's provided
    if inv is not None:
        for network in inv:
            for station in network:
                for channel in station:
                    stats = dict(
                        network=network.code,
                        station=station.code,
                        location=channel.location_code,
                        channel=channel.code,
                        sampling_rate=dummy_sampling_rate,
                    )
                    st += Trace(
                        data=np.zeros(
                            int(data_length_seconds * stats['sampling_rate'])
                        ),
                        header=stats,
                    )
        for tr in st:
            coords = inv.get_coordinates(tr.id)
            tr.stats.longitude = coords['longitude']
            tr.stats.latitude = coords['latitude']

    # Add Traces from fake_station_dict, if it's provided
    if fake_station_dict is not None:
        for station, coords in fake_station_dict.items():
            for channel in 'HXZ', 'HXN', 'HXE':
                stats = dict(
                    network='XX',
                    station=station,
                    channel=channel,
                    location='SE',  # Synthetics Engine
                    sampling_rate=dummy_sampling_rate,
                    longitude=coords[1],
                    latitude=coords[0],
                )
                st += Trace(
                    data=np.zeros(int(data_length_seconds * stats['sampling_rate'])),
                    header=stats,
                )
    if fake_start_time is not None:
        for tr in st:
            tr.stats.starttime = fake_start_time

    # Create LSData instance
    instance = LSData(
        st,
        source_lat=source_lat,
        source_lon=source_lon,
        remove_response=False,
        skip_zne_rotation=skip_zne_rotation,
    )

    return instance


def _rotate_to_rtz(st, skip_zne_rotation):
    r"""Rotate all components of a Stream into radial–transverse–vertical.

    This function performs a two-step rotation. First, it rotates all channels into ZNE
    (vertical, north, west) using irisws-fedcatalog metadata (even stations that claim
    to already be in ZNE). Then, it rotates into RTZ (radial, transverse, vertical).

    Args:
        st (:class:`~obspy.core.stream.Stream`): Input Stream to be rotated with
            ``stats.back_azimuth`` defined
        skip_zne_rotation (bool): If `True`, then the ->ZNE rotation step is skipped.
            This is a necessary flag if the stations used do not have metadata in the
            irisws-fedcatalog (e.g., for synthetic cases)

    Returns:
        :class:`~obspy.core.stream.Stream`: Rotated Stream
    """

    st_rot = st.copy()  # Work on a copy of the data

    # Assemble info for inventory gather
    networks = np.unique([tr.stats.network for tr in st_rot])
    stations = np.unique([tr.stats.station for tr in st_rot])
    channels = np.unique([tr.stats.channel for tr in st_rot])

    # If we're doing the full ->ZNE rotation step prior to NE->RT
    if not skip_zne_rotation:
        # Grab inventory for orientation info
        client = RoutingClient('iris-federator')
        inv = client.get_stations(
            network=','.join(networks),
            station=','.join(stations),
            channel=','.join(channels),
            starttime=st_rot[0].stats.starttime,
            endtime=st_rot[0].stats.endtime,
            level='channel',
        )

        # Rotate to ZNE (must put 'ZNE' in list of components to handle ZNE!)
        st_rot.rotate('->ZNE', components=['ZNE', 'Z12', '123'], inventory=inv)

        # Report on what was modified by the above call
        for tr, tr_rot in zip(st.copy().sort(), st_rot.copy().sort()):
            try:
                np.testing.assert_allclose(tr_rot.data, tr.data, verbose=True)
            except AssertionError as error:
                orientation = inv.get_orientation(tr.id)
                message = '{} -> {}\n{}\n{}\n'.format(
                    tr.id,
                    tr_rot.id,
                    orientation,
                    error.__str__().strip().replace('\n\n', '\n'),
                )
                print(message)

    # Rotate to ZRT on a station-by-station basis
    for station in stations:
        st_sta = st_rot.select(station=station)  # Just select Traces for this station
        st_sta.rotate('NE->RT')

    return st_rot
