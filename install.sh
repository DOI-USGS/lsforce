#!/bin/bash

# USAGE:
#   bash install.sh    # Standard (user) install
#   bash install.sh 1  # Developer install

platform=$(uname)
if [ "$platform" == 'Linux' ]
then
    profile=~/.bashrc
    miniconda_url=https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
elif [ "$platform" == 'FreeBSD' ] || [ "$platform" == 'Darwin' ]
then
    profile=~/.bash_profile
    miniconda_url=https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
else
    echo 'Unsupported platform. Exiting.'
    exit 1
fi

# Name of package
PACKAGE_NAME=lsforce
# Name of new conda environment
ENV_NAME=$PACKAGE_NAME

# Is conda installed?
if ! conda --version
then
    echo 'No conda detected, installing miniconda...'

    # Try to download shell script using curl
    if ! curl -L $miniconda_url -o miniconda.sh
    then
        echo 'Failed to create download miniconda installer shell script. Exiting.'
        exit 1
    fi

    # Try to install miniconda
    echo 'Install directory: ~/miniconda'
    if ! bash miniconda.sh -f -b -p ~/miniconda
    then
        echo 'Failed to run miniconda installer shell script. Exiting.'
        exit 1
    fi

    # Set up conda activate command, see:
    # https://docs.anaconda.com/anaconda/install/silent-mode/#linux-macos
    eval "$(~/miniconda/bin/conda shell.bash hook)"
    conda init

    # Don't automatically activate the environment
    conda config --set auto_activate_base false

else
    echo "conda detected, installing the $ENV_NAME environment..."
fi

# This is needed to ensure that environment activation works
source $profile

# Try to activate the base conda environment
echo 'Activating the base environment'
if ! conda activate base
then
    echo '"conda activate" failed, trying "source activate" instead...'
    if ! source activate base
    then
        echo 'Failed to activate the base environment. Exiting.'
        exit 1
    fi
fi

# Remove existing environment if it exists
conda remove --yes --name $ENV_NAME --all

# Standard package list:
PACKAGE_LIST=(
    'cartopy'
    'obspy'
    'pyqt'
    'scikit-learn'
    'xarray'
)

# Additional developer packages:
DEVELOPER_PACKAGES=(
    'black'
    'ipython'
    'isort'
    'pytest-cov'
    'recommonmark'
    'sphinx'
    'sphinx_rtd_theme'
    'sphinxcontrib-apidoc'
)

# If user supplied the developer flag, add developer packages to package list
# TODO: Remove `|| true` to actually use this flag. Currently we always install as dev!
if [ "$1" == 1 ] || true
then
    PACKAGE_LIST=( "${PACKAGE_LIST[@]}" "${DEVELOPER_PACKAGES[@]}" )
    echo 'Installing developer packages:'
    echo "${DEVELOPER_PACKAGES[@]}"
fi

# Try to create a conda environment
echo "Creating the $ENV_NAME environment"
if ! conda create --yes --name $ENV_NAME --channel conda-forge "${PACKAGE_LIST[@]}"
then
    echo 'Failed to create conda environment. Resolve any conflicts, then try again.'
    exit 1
fi

# Try to activate the new conda environment
echo "Activating the $ENV_NAME environment"
if ! conda activate $ENV_NAME
then
    echo '"conda activate" failed, trying "source activate" instead...'
    if ! source activate $ENV_NAME
    then
        echo "Failed to activate the $ENV_NAME conda environment. Exiting."
        exit 1
    fi
fi

# Try to upgrade pip, mostly so pip doesn't complain about not being new...
if ! pip install --upgrade pip
then
    echo 'Failed to upgrade pip. Trying to continue...'
fi

# Try to install this package
echo
echo "Installing $PACKAGE_NAME"
if ! pip install --editable .
then
    echo 'Failed to pip install this package. Exiting.'
    exit 1
fi

# Tell user to install CPS
echo
echo 'This code requires Computer Programs in Seismology (CPS), available at:'
echo 'http://www.eas.slu.edu/eqc/eqccps.html'
echo 'You will need to add it to your PATH after installing.'
