lsforce
=======

🚨 _lsforce_ is currently under rapid development. Use at your own risk! 🚨

Installation
------------

1. Install [Computer Programs in Seismology](http://www.eas.slu.edu/eqc/eqccps.html),
   and ensure it's on your `PATH`

2. Create a [conda](https://docs.conda.io/en/latest/) environment:
   ```
   conda create -n lsforce -c conda-forge black cartopy ipython obspy scikit-learn xarray
   ```

3. Install [_waveform_collection_](https://github.com/uafgeotools/waveform_collection)
   into this environment:
   ```
   conda activate lsforce
   pip install git+https://github.com/uafgeotools/waveform_collection.git
   ```

4. Clone this repo:
   ```
   git clone https://code.usgs.gov/ghsc/lhp/lsforce.git
   cd lsforce
   ```

5. Set up an arbitrary run directory and grab model file:
   ```
   mkdir meow
   cd meow
   curl -O http://www.eas.slu.edu/eqc/eqc_cps/TUTORIAL/SPHERICITY/AK135/tak135sph.mod
   ```
