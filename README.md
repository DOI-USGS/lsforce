lsforce
=======

🚨 _lsforce_ is currently under rapid development. Use at your own risk! 🚨

Installation
------------

The following has only been tested on macOS Mojave.

1. Install
   [Computer Programs in Seismology (CPS)](http://www.eas.slu.edu/eqc/eqccps.html), and
   ensure it's on your `PATH`:
   
   * Install [GCC](https://gcc.gnu.org/) with e.g. [Homebrew](https://brew.sh/):
     ```
     brew install gcc
     ```
   * Complete the
     [CPS license form](http://www.eas.slu.edu/eqc/eqc_cps/CPS/cpslisc.html), download
     the resulting archive, and unzip
   * Move the directory `PROGRAMS.330/` to where you'd like to install, then:
     ```
     cd PROGRAMS.330
     ./Setup OSX40
     ./C
     ```
   * Add the executables to your `PATH` by adding the following line to your e.g.
     `~/.bash_profile`:
     ```
     export PATH="$PATH:/path/to/PROGRAMS.330/bin"
     ```
     
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

Usage
-----

A usage example is given in `example.py`.
