lsforce
=======

🚨 _lsforce_ is currently under rapid development. Use at your own risk! 🚨

Installation
------------

Grab model file:
```
curl -O http://www.eas.slu.edu/eqc/eqc_cps/TUTORIAL/SPHERICITY/AK135/tak135sph.mod
```

Create conda environment:
```
conda create -n lsforce -c conda-forge black ipython obspy
```

Install [_waveform_collection_](https://github.com/uafgeotools/waveform_collection) into
this environment:
```
conda activate lsforce
pip install git+https://github.com/uafgeotools/waveform_collection.git
```
