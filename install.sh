#!/usr/bin/env bash

# USAGE:
#   bash install.sh    # Standard (user) install
#   bash install.sh 1  # Developer install

platform=$(uname)
if [ "$platform" == 'Linux' ]
then
    miniforge_url="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-$(uname -m).sh"
elif [ "$platform" == 'Darwin' ]
then
    miniforge_url="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-$(uname -m).sh"
else
    echo 'Unsupported platform. Exiting.'
    exit 1
fi

# Name of package
PACKAGE_NAME=lsforce

# Name of new environment
ENV_NAME=$PACKAGE_NAME

# Is mamba installed?
if ! mamba --version
then
    echo 'No mamba detected.'
    # Is conda installed?
    if ! conda --version
    then
        echo 'No conda detected, either — installing Miniforge...'

        # Try to download shell script using curl
        if ! curl -L $miniforge_url -o Miniforge3.sh
        then
            echo 'Failed to download Miniforge installer shell script. Exiting.'
            exit 1
        fi

        # Try to install Miniforge
        echo 'Install directory: ~/miniforge3'
        if ! bash Miniforge3.sh -f -b -p ~/miniforge3
        then
            echo 'Failed to run Miniforge installer shell script. Exiting.'
            exit 1
        fi

        # "create the path to conda" — see:
        # https://github.com/conda-forge/miniforge/blob/main/README.md#as-part-of-a-ci-pipeline
        source ~/miniforge3/etc/profile.d/conda.sh

        echo 'mamba installed.'
        CONDA_CMD=mamba
    else
        echo 'conda detected.'
        CONDA_CMD=conda
    fi
else
    echo 'mamba detected.'
    CONDA_CMD=mamba
fi

# "create the path to conda" — see:
# https://github.com/conda-forge/miniforge/blob/main/README.md#as-part-of-a-ci-pipeline
source $(conda info --base)/etc/profile.d/conda.sh

# Try to activate the base environment
echo 'Activating the base environment.'
if ! conda activate base
then
    echo 'Failed to activate the base environment. Exiting.'
    exit 1
fi

# Remove existing environment if it exists
conda remove --yes --name $ENV_NAME --all

# Try to create the environment
echo "Creating the $ENV_NAME environment"
if ! $CONDA_CMD env create --yes --name $ENV_NAME --file environment.yml
then
    echo 'Failed to create environment. Exiting.'
    exit 1
fi

# Try to activate the new environment
echo "Activating the $ENV_NAME environment."
if ! conda activate $ENV_NAME
then
    echo '"conda activate" failed, trying "source activate" instead...'
    if ! source activate $ENV_NAME
    then
        echo "Failed to activate the $ENV_NAME environment. Exiting."
        exit 1
    fi
fi

# Try to upgrade pip, mostly so pip doesn't complain about not being new...
if ! pip install --upgrade pip
then
    echo 'Failed to upgrade pip. Trying to continue...'
fi

# If user supplied the developer flag, install developer packages
if [ "$1" == 1 ]
then
    echo 'Installing developer packages...'
    if ! pip install --requirement requirements.txt
    then
        echo 'Failed to install developer packages. Exiting.'
    fi
fi

# Try to install this package
echo
echo "Installing $PACKAGE_NAME using pip."
if ! pip install --editable .
then
    echo "Failed to pip install $PACKAGE_NAME. Exiting."
    exit 1
fi
