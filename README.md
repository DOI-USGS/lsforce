lsforce
=======

[![pipeline status](https://code.usgs.gov/ghsc/users/ltoney/lsforce/badges/master/pipeline.svg)](https://code.usgs.gov/ghsc/users/ltoney/lsforce/pipelines/latest)
[![coverage report](https://code.usgs.gov/ghsc/users/ltoney/lsforce/badges/master/coverage.svg)](https://code.usgs.gov/ghsc/users/ltoney/lsforce/-/jobs)

🚨 _lsforce_ is currently under rapid development. Use at your own risk! 🚨

Installation
------------

The following has only been tested on macOS Mojave.

1. Install
   [Computer Programs in Seismology (CPS)](http://www.eas.slu.edu/eqc/eqccps.html), and
   ensure it's on your `PATH`:

   * Install [GCC](https://gcc.gnu.org/) with e.g. [Homebrew](https://brew.sh/):
     ```shell
     brew install gcc
     ```
   * Complete the
     [CPS license form](http://www.eas.slu.edu/eqc/eqc_cps/CPS/cpslisc.html), download
     the resulting archive, and unzip
   * Move the directory `PROGRAMS.330` to where you'd like to install, then:
     ```shell
     cd PROGRAMS.330
     ./Setup OSX40
     ./C
     ```
   * Add the executables to your `PATH` by adding the following line to your e.g.
     `~/.bash_profile`:
     ```shell
     export PATH="$PATH:/path/to/PROGRAMS.330/bin"
     ```

2. Clone this repo and run the installation script, which creates a
   [conda](https://docs.conda.io/en/latest/) environment named `lsforce` and installs
   the _lsforce_ package into the environment:
   ```shell
   git clone https://code.usgs.gov/ghsc/lhp/lsforce.git
   cd lsforce
   bash install.sh  # Or `bash install.sh 1` if you want developer tools as well
   ```

4. Set up an arbitrary run directory and grab model file:
   ```shell
   mkdir meow
   cd meow
   curl -O http://www.eas.slu.edu/eqc/eqc_cps/TUTORIAL/SPHERICITY/AK135/tak135sph.mod
   ```

Documentation
-------------

A usage example is given in `example.py`.

To build the documentation, first ensure that you installed the developer tools (`bash
install.sh 1`), which are required for documentation building. Then:
```shell
conda activate lsforce
cd doc
make html
open _build/html/index.html  # macOS command to open file in browser
```

Testing
-------

Tests are located in the `tests` directory. To run the tests, first ensure that you
installed the developer tools (`bash install.sh 1`), which are required for testing.
Then:
```shell
conda activate lsforce
pytest --cov=lsforce tests/
```
