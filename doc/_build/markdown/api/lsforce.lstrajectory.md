# lsforce.lstrajectory module


### _class_ lsforce.lstrajectory.LSTrajectory(force, mass=None, target_length=None, duration=None, detrend_velocity=None, zeroacc=None)
Bases: [`object`](https://docs.python.org/3/library/functions.html#object)

Class for force inversion derived trajectories.


#### force()
Inversion results used to compute
this trajectory


* **Type**

    [`LSForce`](lsforce.lsforce.md#lsforce.lsforce.LSForce)



#### mass_requested()
[kg] Mass specified


* **Type**

    [int](https://docs.python.org/3/library/functions.html#int) or [float](https://docs.python.org/3/library/functions.html#float)



#### mass_actual()
[kg] Mass used (same as mass_requested if target_length
is not specified)


* **Type**

    [int](https://docs.python.org/3/library/functions.html#int)



#### target_length()
[km] Center-of-mass runout length of event, None
if not specified


* **Type**

    [int](https://docs.python.org/3/library/functions.html#int) or [float](https://docs.python.org/3/library/functions.html#float)



#### jackknife()
Jackknifed
trajectory results


* **Type**

    [`AttribDict`](https://docs.obspy.org/packages/autogen/obspy.core.util.attribdict.AttribDict.html#obspy.core.util.attribdict.AttribDict)



#### acceleration()
[m^2/s] Computed
acceleration with Z, E, N components as attributes


* **Type**

    [`AttribDict`](https://docs.obspy.org/packages/autogen/obspy.core.util.attribdict.AttribDict.html#obspy.core.util.attribdict.AttribDict)



#### velocity()
[m/s] Computed
velocity with Z, E, N components as attributes


* **Type**

    [`AttribDict`](https://docs.obspy.org/packages/autogen/obspy.core.util.attribdict.AttribDict.html#obspy.core.util.attribdict.AttribDict)



#### displacement()
[m] Computed
displacement with Z, E, N components as attributes


* **Type**

    [`AttribDict`](https://docs.obspy.org/packages/autogen/obspy.core.util.attribdict.AttribDict.html#obspy.core.util.attribdict.AttribDict)



#### horizontal_distance()
[m] Computed horizontal distance


* **Type**

    [`ndarray`](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)



#### traj_tvec()
[s] Time array for all trajectory arrays


* **Type**

    [`ndarray`](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray)


Create an LSTrajectory object.


* **Parameters**

    
    * **force** ([`LSForce`](lsforce.lsforce.md#lsforce.lsforce.LSForce)) – Completed force inversion


    * **mass** ([*int*](https://docs.python.org/3/library/functions.html#int)* or *[*float*](https://docs.python.org/3/library/functions.html#float)) – [kg] Mass of event. If None, the mass is computed
    using target_length, which must be specified


    * **target_length** ([*int*](https://docs.python.org/3/library/functions.html#int)* or *[*float*](https://docs.python.org/3/library/functions.html#float)) – [km] Center-of-mass runout length of event. If
    None, mass must be specified


    * **duration** ([*int*](https://docs.python.org/3/library/functions.html#int)* or *[*float*](https://docs.python.org/3/library/functions.html#float)) – [s] If not None, only use the force time series
    from 0–duration seconds in the trajectory calculation


    * **detrend_velocity** – [s] If provided, require the velocity to linearly go to
    zero at this time; if None, don’t detrend


    * **zeroacc** – [s] If provided, require the acceleration to be zero after
    this time, usually when forces are indistinguishable from
    zero, to reduce noise at end of trajectory. This is best used
    with the detrend_velocity option to avoid allowing the landslide to
    slide forever.



#### plot_trajectory(elevation_profile=False, plot_jackknife=False, image=None, dem=None, reference_point=None)
Plot trajectory results with context.


* **Parameters**

    
    * **elevation_profile** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – If True, plot vertical displacement versus
    horizontal runout distance ($H$ vs. $L$) instead of a map
    view


    * **plot_jackknife** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – Toggle plotting jackknifed displacements as well (if
    available)


    * **image** ([`DataArray`](https://docs.xarray.dev/en/stable/generated/xarray.DataArray.html#xarray.DataArray)) – An image with coordinates defined in km
    with the origin (0, 0) being the start location of the trajectory


    * **dem** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – A UTM-projected DEM GeoTIFF to slice thru for elevation profile
    plot


    * **reference_point** ([*int*](https://docs.python.org/3/library/functions.html#int)* or *[*float*](https://docs.python.org/3/library/functions.html#float)* or *[*list*](https://docs.python.org/3/library/stdtypes.html#list)) – If not None, plot a dot on
    trajectory, and line on colorbar, at this specified time(s) for
    reference



* **Returns**

    Output figure handle



* **Return type**

    [`Figure`](https://matplotlib.org/stable/api/figure_api.html#matplotlib.figure.Figure)
