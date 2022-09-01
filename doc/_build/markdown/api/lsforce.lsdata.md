# lsforce.lsdata module


### _class_ lsforce.lsdata.LSData(st, source_lat, source_lon, remove_response=True, skip_zne_rotation=False)
Bases: [`object`](https://docs.python.org/3/library/functions.html#object)

Class for force inversion data that is an extension of an ObsPy Stream.


#### st_orig()
Original input Stream st.


* **Type**

    [`Stream`](https://docs.obspy.org/packages/autogen/obspy.core.stream.Stream.html#obspy.core.stream.Stream)



#### st_proc()
Stream rotated into RTZ (radial,
transverse, vertical) relative to source_lat, source_lon.


* **Type**

    [`Stream`](https://docs.obspy.org/packages/autogen/obspy.core.stream.Stream.html#obspy.core.stream.Stream)



#### source_lat()
Latitude in decimal degrees of centroid of landslide
source location.


* **Type**

    [float](https://docs.python.org/3/library/functions.html#float)



#### source_lon()
Longitude in decimal degrees of centroid of landslide
source location.


* **Type**

    [float](https://docs.python.org/3/library/functions.html#float)


Create an LSData object.


* **Parameters**

    
    * **st** ([`Stream`](https://docs.obspy.org/packages/autogen/obspy.core.stream.Stream.html#obspy.core.stream.Stream)) – Stream object with
    `tr.stats.latitude` and `tr.stats.longitude` defined and station
    response info attached to each trace in the Stream.


    * **source_lat** ([*float*](https://docs.python.org/3/library/functions.html#float)) – Latitude in decimal degrees of centroid of landslide
    source location


    * **source_lon** ([*float*](https://docs.python.org/3/library/functions.html#float)) – Longitude in decimal degrees of centroid of landslide
    source location


    * **remove_response** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – Correct for station response to displacement units.
    Set to False to handle response removal manually at an earlier step.


    * **skip_zne_rotation** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – If True, then the ->ZNE rotation step is
    skipped. This is a necessary flag if the stations used do not have
    metadata in the irisws-fedcatalog (e.g., for synthetic cases)



#### plot_data(equal_scale=True, period_range=None)
Create a record section plot of waveforms in st_proc.


* **Parameters**

    
    * **equal_scale** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – If True, all plots will share the same y-axis scale


    * **period_range** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)* or *[*tuple*](https://docs.python.org/3/library/stdtypes.html#tuple)) – If not None, filter the data between
    period_range[0] and period_range[1], given in seconds



* **Returns**

    Output figure handle



* **Return type**

    [`Figure`](https://matplotlib.org/stable/api/figure_api.html#matplotlib.figure.Figure)



#### plot_stations(region=None, label_stations=False, gshhs_scale='auto')
Create a map showing stations and event location.


* **Parameters**

    
    * **region** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)* or *[*tuple*](https://docs.python.org/3/library/stdtypes.html#tuple)) – Array of the form [lonmin, lonmax, latmin, latmax]
    specifying the desired map region in decimal degrees. If None, we
    automatically pick a region that includes the event and stations


    * **label_stations** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – If True, label stations with their codes


    * **gshhs_scale** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – Resolution for coastlines; one of ‘auto’, ‘coarse’,
    ‘low’, ‘intermediate’, ‘high’, or ‘full’



* **Returns**

    Output figure handle



* **Return type**

    [`Figure`](https://matplotlib.org/stable/api/figure_api.html#matplotlib.figure.Figure)



### lsforce.lsdata.make_lsdata_syn(inv, fake_station_dict, source_lat, source_lon, data_length_seconds)
Wrapper which creates an `LSData` object for
forward modeling applications.


* **Parameters**

    
    * **inv** ([`Inventory`](https://docs.obspy.org/packages/autogen/obspy.core.inventory.inventory.Inventory.html#obspy.core.inventory.inventory.Inventory)) – ObsPy Inventory
    object containing the **real** stations for which synthetics should be
    computed


    * **fake_station_dict** ([*dict*](https://docs.python.org/3/library/stdtypes.html#dict)) – Dictionary with keys specifying station names, and
    values as two-element lists [latitude, longitude], of **fake** stations
    for which synthetics should be computed


    * **source_lat** ([*float*](https://docs.python.org/3/library/functions.html#float)) – Latitude in decimal degrees of centroid of landslide
    source location


    * **source_lon** ([*float*](https://docs.python.org/3/library/functions.html#float)) – Longitude in decimal degrees of centroid of landslide
    source location


    * **data_length_seconds** ([*int*](https://docs.python.org/3/library/functions.html#int)* or *[*float*](https://docs.python.org/3/library/functions.html#float)) – [s] Length of synthetic data
