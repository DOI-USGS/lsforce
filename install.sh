#!/bin/bash

unamestr=`uname`
if [ "$unamestr" == 'Linux' ]; then
    prof=~/.bashrc
    mini_conda_url=https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
elif [ "$unamestr" == 'FreeBSD' ] || [ "$unamestr" == 'Darwin' ]; then
    prof=~/.bash_profile
    mini_conda_url=https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
else
    echo "Unsupported environment. Exiting."
    exit
fi

source $prof

# Name of package
PACKAGE_NAME=lsforce
# Name of new environment
ENV_NAME=$PACKAGE_NAME
# Python version
PYTHON_VERSION=3.8

# Set to 1 if you are a developer and want Black, IPython, etc. installed
DEVELOPER=0

# Is conda installed?
conda --version
if [ $? -ne 0 ]; then
    echo "No conda detected, installing miniconda..."

    curl -L $mini_conda_url -o miniconda.sh;

    # If curl fails, bow out gracefully
    if [ $? -ne 0 ];then
        echo "Failed to create download miniconda installer shell script. Exiting."
        exit 1
    fi

    echo "Install directory: $HOME/miniconda"

    bash miniconda.sh -f -b -p $HOME/miniconda

    # If miniconda.sh fails, bow out gracefully
    if [ $? -ne 0 ];then
        echo "Failed to run miniconda installer shell script. Exiting."
        exit 1
    fi

    . $HOME/miniconda/etc/profile.d/conda.sh
else
    echo "conda detected, installing $ENV_NAME environment..."
fi

# Add source command to profile file if it isn't already there
grep -q "/etc/profile.d/conda.sh" $prof
if [ $? -ne 0 ]; then
    echo ". $_CONDA_ROOT/etc/profile.d/conda.sh" >> $prof
fi

# Start in conda base environment
echo "Activating base environment"
conda activate base

# Remove existing environment if it exists
conda remove --yes --name $ENV_NAME --all

dev_list=(
    "black"
    "ipython"
)

# Package list:
package_list=(
    "python=$PYTHON_VERSION"
    "cartopy"
    "pyqt"
    "scikit-learn"
    "xarray"
)

if [ $DEVELOPER == 1 ]; then
    package_list=( "${package_list[@]}" "${dev_list[@]}" )
    echo ${package_list[*]}
fi

# Create a conda environment
echo "Creating the $ENV_NAME environment"
conda create --yes --name $ENV_NAME --channel conda-forge ${package_list[*]}

# Bail out at this point if the conda create command fails.
# Clean up zip files we've downloaded
if [ $? -ne 0 ]; then
    echo "Failed to create conda environment. Resolve any conflicts, then try again."
    exit
fi

# Activate the new environment
echo "Activating the $ENV_NAME environment"
conda activate $ENV_NAME

# If conda activate fails, bow out gracefully
if [ $? -ne 0 ];then
    echo "Failed to activate $ENV_NAME conda environment. Exiting."
    exit 1
fi

# Upgrade pip, mostly so pip doesn't complain about not being new...
pip install --upgrade pip

# If pip upgrade fails, complain but try to keep going
if [ $? -ne 0 ];then
    echo "Failed to upgrade pip, trying to continue..."
    exit 1
fi

# Install this package
echo "Installing $PACKAGE_NAME"
pip install --editable .

# If pip install fails, bow out gracefully
if [ $? -ne 0 ];then
    echo "Failed to pip install this package. Exiting."
    exit 1
fi

# Tell user to install CPS
echo "This code requires Computer Programs in Seismology (CPS), available at:"
echo "http://www.eas.slu.edu/eqc/eqccps.html"
echo "You'll need to add it to your PATH after installing."
