from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import numpy as np


def stationmap_coords(lats, lons, stations=None, channels=None, event_lat=None, event_lon=None, chanshow=''):
    """Plot stations specified without having to specify a source location or pull info from database

    Args:
        lats: np array of latitudes corresponding to stations
        lons: np array of longitudes
        stations: list of station names
        channels: list of channel names corresponding to each station names
        event_lat: latitude of event to plot relative to
        event_lon: ditto
        chanshow - channels to plot, '' or '*' plots all stations specified, else use form chanshow = 'ENZ,EHZ,HHZ'

    Returns:
        figure handle
    """
    fig = plt.figure()
    lats = np.array(lats)
    lons = np.array(lons)
    if event_lat is None:
        event_lat = np.mean(lats)
        event_lon = np.mean(lons)
    alllats = np.hstack([lats, event_lat])
    alllons = np.hstack([lons, event_lon])
    difflat = max(alllats)-min(alllats)
    difflon = max(alllons)-min(alllons)
    m = Basemap(projection='stere', lat_0=event_lat, lon_0=event_lon, llcrnrlat=min(alllats)-0.1*difflat, urcrnrlat=max(alllats)+0.1*difflat, llcrnrlon=min(alllons)-0.1*difflon, urcrnrlon=max(alllons)+0.1*difflon, resolution='i')
    m.shadedrelief(alpha=0.4)
    m.drawcoastlines(color='#6E6E6E')
    m.drawstates(color='#6E6E6E')
    m.drawcountries(color='#6E6E6E')
    m.drawparallels(np.arange(np.floor(min(alllats)), np.ceil(max(alllats)), 1), labels=[1, 1, 0, 0])
    m.drawmeridians(np.arange(np.floor(min(alllons)), np.ceil(max(alllons)), 1), labels=[0, 0, 0, 1])
    #m.fillcontinents(color='white',lake_color='#ccFFFF')
    m.drawmapboundary(fill_color='#ccFFFF')
    #m.scatter(lons,lats,latlon=True)
    xs, ys = m(lons, lats)
    xyplotted = []
    yoffset = 0.04*(m.ymax-m.ymin)
    dmin = yoffset
    if event_lat is not None:
        m.scatter(event_lon, event_lat, c='g', s=800, latlon=True, marker='*', zorder=9)
    v = 0
    BB = []
    INF = []
    SP = []
    UK = []
    if stations is not None and channels is not None:
        for x, y, sta, chan in zip(xs, ys, stations, channels):
            if chan in chanshow or chanshow is '' or chanshow is '*':
                if 'BH' in chan or 'HH' in chan:
                    BB = plt.scatter(x, y, c='r', marker='^', s=75, zorder=10+v)
                    v += 1
                elif 'DF' in chan:
                    INF = plt.scatter(x, y, c='m', marker='s', s=75, zorder=10+v)
                    v += 1
                else:
                    SP = plt.scatter(x, y, marker='o', s=50, zorder=10+v)
                    v += 1
                dist = [np.sqrt((x-x0)**2+(y-y0)**2) for x0, y0 in xyplotted]
                if not dist or min(dist) > dmin:
                    plt.text(x, y, sta, ha='right', va='top')
                xyplotted.append((x, y))
    else:
        for x, y in zip(xs, ys):
            UK = plt.scatter(x, y, c='k', marker='^', s=75, zorder=10+v)
            v += 1
            dist = [np.sqrt((x-x0)**2+(y-y0)**2) for x0, y0 in xyplotted]
            if not dist or min(dist) > dmin:
                plt.text(x, y, sta, ha='right', va='top')
                xyplotted.append((x, y))
    plt.legend([SP, BB, INF, UK], ['Short period', 'Broadband', 'Infrasound', 'Unknown'])
    plt.show()
    return fig
