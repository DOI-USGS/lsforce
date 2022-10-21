# lsforce.lsforce module


### _class_ lsforce.lsforce.LSForce(data, data_sampling_rate, main_folder=None, method='full')
Bases: [`object`](https://docs.python.org/3/library/functions.html#object)

Class for performing force inversions.


#### gf_dir()
Directory containing Green’s functions


* **Type**

    [str](https://docs.python.org/3/library/stdtypes.html#str)



#### gf_computed()
Whether or not Green’s functions have been computed for this
object


* **Type**

    [bool](https://docs.python.org/3/library/functions.html#bool)



#### filtered_gf_st()
Stream containing filtered
Green’s functions


* **Type**

    [`Stream`](https://docs.obspy.org/packages/autogen/obspy.core.stream.Stream.html#obspy.core.stream.Stream)



#### inversion_complete()
Whether or not the inversion has been run


* **Type**

    [bool](https://docs.python.org/3/library/functions.html#bool)



#### filter()
Dictionary with keys `'freqmin'`, `'freqmax'`,
`'zerophase'`, `'periodmin'`, `'periodmax'`, and `'order'`
specifying filter parameters


* **Type**

    [dict](https://docs.python.org/3/library/stdtypes.html#dict)



#### data_length()
Length in samples of each data trace


* **Type**

    [int](https://docs.python.org/3/library/functions.html#int)



#### force_sampling_rate()
[Hz] The sampling rate of the force-time
function


* **Type**

    [int](https://docs.python.org/3/library/functions.html#int) or [float](https://docs.python.org/3/library/functions.html#float)



#### W()
Weight matrix


* **Type**

    2D array



#### Wvec()
Weight vector


* **Type**

    1D array



#### jackknife()
Dictionary with
keys `'Z'`, `'N'`, `'E'`, `'VR_all'`, `'alphas'`, `'num_iter'`,
and `'frac_delete'` containing jackknife parameters and results


* **Type**

    [`AttribDict`](https://docs.obspy.org/packages/autogen/obspy.core.util.attribdict.AttribDict.html#obspy.core.util.attribdict.AttribDict)



#### angle_magnitude()
Dictionary
with keys `'magnitude'`, `'magnitude_upper'`, `'magnitude_lower'`,
`'vertical_angle'`, and `'horizontal_angle'` containing inversion angle
and magnitude information


* **Type**

    [`AttribDict`](https://docs.obspy.org/packages/autogen/obspy.core.util.attribdict.AttribDict.html#obspy.core.util.attribdict.AttribDict)



#### G()
Design matrix


* **Type**

    2D array



#### d()
Data vector


* **Type**

    1D array



#### model()
Model vector of concatenated components (n x 1) of solution


#### Z()
[N] Vertical force time series extracted from model (positive up)


#### N()
[N] North force time series extracted from model (positive north)


#### E()
[N] East force time series extracted from model (positive east)


#### tvec()
[s] Time vector for forces, referenced using zero_time (if specified)


#### VR()
[%] Variance reduction. Rule of thumb: This should be ~50–80%, if ~100%,
solution is fitting data exactly and results are suspect. If ~5%, model may
be wrong or something else may be wrong with setup


#### dtorig()
Original data vector


#### dtnew()
Modeled data vector (Gm-d)


#### alpha()
Regularization parameter that was used


#### alphafit()
Dictionary with keys `'alphas'`, `'fit'`, and `'size'`
specifying regularization parameters tested


* **Type**

    [dict](https://docs.python.org/3/library/stdtypes.html#dict)


Create an LSForce object.


* **Parameters**

    
    * **data** ([`LSData`](lsforce.lsdata.md#lsforce.lsdata.LSData)) – LSData object, corrected for station
    response but not filtered


    * **data_sampling_rate** ([*int*](https://docs.python.org/3/library/functions.html#int)* or *[*float*](https://docs.python.org/3/library/functions.html#float)) – [Hz] Samples per second to use in
    inversion. All data will be resampled to this rate, and Green’s
    functions will be created with this rate


    * **main_folder** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – If None, will use current folder


    * **method** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – How to parameterize the force-time function. One of ‘full’
    — full inversion using Tikhonov regularization (L2 norm minimization) or
    ‘triangle’ — inversion parameterized using overlapping triangles,
    variation on method of Ekström & Stark (2013)



#### forward(Z, N, E)
Execute the forward problem $\mathbf{d} = \mathbf{G}\mathbf{m}$
using a user-supplied force time series $\mathbf{m}$ composed of
components Z, N, and E.


* **Parameters**

    
    * **Z** – [N] Vertical force time series (positive up)


    * **N** – [N] North force time series (positive north)


    * **E** – [N] East force time series (positive east)



* **Returns**

    [m] Stream containing synthetic
    data, $\mathbf{d}$



* **Return type**

    [`Stream`](https://docs.obspy.org/packages/autogen/obspy.core.stream.Stream.html#obspy.core.stream.Stream)



#### invert(zero_time=None, impose_zero_start=False, add_to_zero=False, duration=None, jackknife=False, num_iter=200, frac_delete=0.5, alpha=None, zero_scaler=2.0, zero_start_taper_length=0, tikhonov_ratios=(1.0, 0.0, 0.0), jk_refine_alpha=False, save_matrices=False)
Performs single-force inversion using Tikhonov regularization.


* **Parameters**

    
    * **zero_time** ([*int*](https://docs.python.org/3/library/functions.html#int)* or *[*float*](https://docs.python.org/3/library/functions.html#float)) – [s] Optional estimated start time of real
    (landslide-related) part of signal, in seconds from start time of
    seismic data. Useful for making figures showing selected start time and
    also for the impose_zero_start option


    * **impose_zero_start** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – Adds weighting matrix to suggest that forces tend
    towards zero prior to zero_time (zero_time must be defined)


    * **add_to_zero** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – Adds weighting matrix to suggest that all components of
    force integrate to zero


    * **duration** ([*int*](https://docs.python.org/3/library/functions.html#int)* or *[*float*](https://docs.python.org/3/library/functions.html#float)) – Maximum duration allowed for the event, starting at
    zero_time if defined, otherwise starting from the beginning of the
    seismic data. Forces after this will tend towards zero. This helps tamp
    down artifacts due to edge effects, etc.


    * **jackknife** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – If True, perform num_iter additional iterations of the
    model while randomly discarding frac_delete of the data


    * **num_iter** ([*int*](https://docs.python.org/3/library/functions.html#int)) – Number of jackknife iterations to perform


    * **frac_delete** ([*int*](https://docs.python.org/3/library/functions.html#int)* or *[*float*](https://docs.python.org/3/library/functions.html#float)) – Fraction (out of 1) of data to discard for each
    iteration, if frac_delete=1, will do leave a one out error analysis


    * **alpha** ([*int*](https://docs.python.org/3/library/functions.html#int)* or *[*float*](https://docs.python.org/3/library/functions.html#float)) – Set regularization parameter. If None, will search
    for best alpha using the L-curve method


    * **zero_scaler** ([*int*](https://docs.python.org/3/library/functions.html#int)* or *[*float*](https://docs.python.org/3/library/functions.html#float)) – Relative strength of zero constraint for
    impose_zero_start and duration options. Ranges from 0 to 10. The
    lower the number, the weaker the constraint. Values up to 30 are
    technically allowed but discouraged because high zero_scaler values
    risk the addition of high frequency oscillations due to the sudden
    release of the constraint


    * **zero_start_taper_length** ([*int*](https://docs.python.org/3/library/functions.html#int)* or *[*float*](https://docs.python.org/3/library/functions.html#float)) – [s] Length of taper for
    impose_zero_start option


    * **tikhonov_ratios** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)* or *[*tuple*](https://docs.python.org/3/library/stdtypes.html#tuple)) – Proportion each regularization method
    contributes to the overall regularization effect, where values
    correspond to [0th order, 1st order, 2nd order]. Must sum to 1


    * **jk_refine_alpha** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – Refine the alpha parameter used for each jackknife
    iteration by searching over order of magnitude around the best alpha
    for the full solution. If False, each jackknife iteration will use the
    same alpha as the main solution (note that this is much faster but can
    result in some jackknife iterations having depressed amplitudes)


    * **save_matrices** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – If True, will save the inverted matrices as
    part of the object (Ghat, dhat, I, L1, L2) in case user wants
    to do additional alpha searching



#### plot_angle_magnitude(xlim=None, ylim=None)
Plot angles and magnitudes of inversion result.


* **Parameters**

    
    * **xlim** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)* or *[*tuple*](https://docs.python.org/3/library/stdtypes.html#tuple)) – x-axis limits


    * **ylim** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)* or *[*tuple*](https://docs.python.org/3/library/stdtypes.html#tuple)) – y-axis limits



* **Returns**

    Output figure handle



* **Return type**

    [`Figure`](https://matplotlib.org/stable/api/figure_api.html#matplotlib.figure.Figure)



#### plot_fits(equal_scale=True, xlim=None)
Create a plot showing the model-produced waveform fit to the data.


* **Parameters**

    
    * **equal_scale** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – If True, all plots will share the same y-axis scale


    * **xlim** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)* or *[*tuple*](https://docs.python.org/3/library/stdtypes.html#tuple)) – [s] Array (length two) of x-axis limits (time relative
    to zero time)



* **Returns**

    Output figure handle



* **Return type**

    [`Figure`](https://matplotlib.org/stable/api/figure_api.html#matplotlib.figure.Figure)



#### plot_forces(subplots=False, xlim=None, ylim=None, same_y=True, highf_tr=None, hfshift=0.0, hfylabel=None, infra_tr=None, infra_shift=0, jackshowall=False)
Plot inversion result.


* **Parameters**

    
    * **subplots** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – If True, make subplots for components, otherwise plot all
    on one plot


    * **xlim** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)* or *[*tuple*](https://docs.python.org/3/library/stdtypes.html#tuple)) – x-axis limits


    * **ylim** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)* or *[*tuple*](https://docs.python.org/3/library/stdtypes.html#tuple)) – y-axis limits


    * **same_y** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – If True, use same y-axis limits for all plots


    * **highf_tr** ([`Trace`](https://docs.obspy.org/packages/autogen/obspy.core.trace.Trace.html#obspy.core.trace.Trace)) – Seismic trace with start time
    identical to start time of the data used in the inversion


    * **hfshift** ([*int*](https://docs.python.org/3/library/functions.html#int)* or *[*float*](https://docs.python.org/3/library/functions.html#float)) – [s] Time shift for seismic trace


    * **hfylabel** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – Label used for seismic trace. If not defined, will use
    station name


    * **infra_tr** ([`Trace`](https://docs.obspy.org/packages/autogen/obspy.core.trace.Trace.html#obspy.core.trace.Trace)) – Infrasound trace with start
    time identical to start time of the data used in the inversion


    * **infra_shift** ([*int*](https://docs.python.org/3/library/functions.html#int)* or *[*float*](https://docs.python.org/3/library/functions.html#float)) – [s] Time shift for infrasound trace


    * **jackshowall** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – If True and jackknife was run, will show all
    individual runs (changes subplots to True)



* **Returns**

    Output figure handle



* **Return type**

    [`Figure`](https://matplotlib.org/stable/api/figure_api.html#matplotlib.figure.Figure)



#### saverun(prefix, filepath=None, timestamp=False, figs2save=None, figs2save_names=None, light=True, filetype='png')
Save a force inversion run for later use.

**WARNING**: Do not expect this to work if you have the `autoreload`
[IPython extension](https://ipython.readthedocs.io/en/stable/config/extensions/autoreload.html) enabled!


* **Parameters**

    
    * **prefix** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – Run name to prepend to all output files


    * **filepath** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – Full path to directory where all files should be saved. If
    None, will use self.main_folder


    * **timestamp** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – Name results with current time to avoid overwriting
    previous results


    * **figs2save** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)* or *[*tuple*](https://docs.python.org/3/library/stdtypes.html#tuple)) – Figure handles to save


    * **figs2save_names** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)* or *[*tuple*](https://docs.python.org/3/library/stdtypes.html#tuple)) – Names of figures (appends to end)


    * **light** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – If True, does not save seismic data with object to save size


    * **filetype** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – Filetype given as extension, e.g. ‘png’



#### setup(period_range, syngine_model=None, cps_model=None, triangle_half_width=None, source_depth=0, weights=None, noise_window_dur=None, filter_order=2, zerophase=True, skip_datafilter=False)
Downloads/computes Green’s functions (GFs) and creates all matrices.


* **Parameters**

    
    * **period_range** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)* or *[*tuple*](https://docs.python.org/3/library/stdtypes.html#tuple)) – [s] Bandpass filter corners


    * **syngine_model** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – Name of Syngine model to use. If this is not None, then
    we calculate GFs using Syngine (preferred)


    * **cps_model** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – Filename of CPS model to use. If this is not None, then we
    calculate GFs using CPS


    * **triangle_half_width** ([*int*](https://docs.python.org/3/library/functions.html#int)* or *[*float*](https://docs.python.org/3/library/functions.html#float)) – [s] Half-width of triangles; only used
    if the triangle method is being used


    * **source_depth** ([*int*](https://docs.python.org/3/library/functions.html#int)* or *[*float*](https://docs.python.org/3/library/functions.html#float)) – [m] Source depth in meters


    * **weights** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)* or *[*tuple*](https://docs.python.org/3/library/stdtypes.html#tuple)* or *[*str*](https://docs.python.org/3/library/stdtypes.html#str)) – If None, no weighting is applied. An array
    of floats with length `st_proc.count()` (and in the order of the `st_proc`
    attribute of the [`LSData`](lsforce.lsdata.md#lsforce.lsdata.LSData) object) applies manual
    weighting. If ‘prenoise’, uses standard deviation of a noise window
    defined by noise_window_dur to weight. If ‘distance’, weights by 1 /
    distance


    * **noise_window_dur** ([*int*](https://docs.python.org/3/library/functions.html#int)* or *[*float*](https://docs.python.org/3/library/functions.html#float)) – [s] Length of noise window for ‘prenoise’
    weighting scheme (if not None, weights is set to ‘prenoise’)


    * **filter_order** ([*int*](https://docs.python.org/3/library/functions.html#int)) – Order of filter applied over period_range


    * **zerophase** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – If True, zero-phase filtering will be used


    * **skip_datafilter** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – If True, filtering will not be applied to
    the input data and will only be applied to the Green’s functions.
    This should be chosen only if the data were pre-filtered
    manually already with the same band as period_range so the
    user doesn’t want to filter them again



#### write_forces(prefix, filepath=None, timestamp=False)
Save E, N, Z forces to a text file for non-*lsforce* users.

File can be read in using, e.g., NumPy as follows:

```python
import numpy as np
e, n, z = np.loadtxt('/path/to/file.txt', unpack=True)
```


* **Parameters**

    
    * **prefix** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – Run name to prepend to file


    * **filepath** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – Full path to directory where file should be saved.
    If None, will use self.main_folder


    * **timestamp** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – Name file with current time to avoid overwriting
    previous results



### lsforce.lsforce.Lcurve(fit1, size1, alphas, bestalpha=None)
Plot an L-curve.


* **Parameters**

    
    * **fit1** (*1D array*) – List of residuals


    * **size1** (*1D array*) – List of model norms


    * **alphas** (*1D array*) – List of alphas tried


    * **bestalpha** ([*float*](https://docs.python.org/3/library/functions.html#float)) – The alpha value chosen



* **Returns**

    figure handle



### lsforce.lsforce.find_alpha(Ghat, dhat, I, L1=0, L2=0, tikhonov_ratios=(1.0, 0.0, 0.0), rough=False, range_rough=None, int_rough=0.75, plot_Lcurve=True)
Finds best regularization (trade-off) parameter alpha.

Computes model with many values of alpha, plots L-curve, and finds point
of steepest curvature where slope is negative.


* **Parameters**

    
    * **Ghat** (*array*) – (m x n) matrix


    * **dhat** (*array*) – (1 x n) array of weighted data


    * **I** (*array*) – Identity matrix


    * **L1** (*array*) – First order roughening matrix. If 0, will use only
    0th-order Tikhonov regularization


    * **L2** (*array*) – Second order roughening matrix. If 0, will use only
    0th-order Tikhonov regularization


    * **tikhonov_ratios** ([*list*](https://docs.python.org/3/library/stdtypes.html#list)* or *[*tuple*](https://docs.python.org/3/library/stdtypes.html#tuple)) – Proportion each regularization method
    contributes to the overall regularization effect, where values
    correspond to [0th order, 1st order, 2nd order]. Must sum to 1


    * **rough** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – If False, will do two iterations to fine tune the alpha
    parameter. The second iteration searches over +/- one order of
    magnitude from the best alpha found from the first round. If
    True, time will be saved because it will only do one round of
    searching.


    * **range_rough** ([*tuple*](https://docs.python.org/3/library/stdtypes.html#tuple)) – Lower and upper bound of range to search over in
    log units. If None, the program will choose a range based on the
    norm of `Ghat`


    * **int_rough** ([*float*](https://docs.python.org/3/library/functions.html#float)) – Interval, in log units, to use for rough alpha search


    * **plot_Lcurve** ([*bool*](https://docs.python.org/3/library/functions.html#bool)) – Toggle showing the L-curve plot



* **Returns**

    Tuple containing:


    * **bestalpha** (float) – The optimal alpha


    * **fit1** (1D array) – List of residuals


    * **size1** (1D array) – List of model norms


    * **alphas** (1D array) – List of alphas tried




* **Return type**

    [tuple](https://docs.python.org/3/library/stdtypes.html#tuple)



### lsforce.lsforce.readrun(filename)
Read in a saved LSForce object.

**WARNING**: Do not expect this to work if you have the `autoreload`
[IPython extension](https://ipython.readthedocs.io/en/stable/config/extensions/autoreload.html) enabled!


* **Parameters**

    **filename** ([*str*](https://docs.python.org/3/library/stdtypes.html#str)) – File path to LSForce object saved using
    `saverun()`



* **Returns**

    Saved LSForce object



* **Return type**

    `LSForce`
