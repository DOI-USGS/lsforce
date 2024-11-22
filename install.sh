#!/bin/bash

# USAGE:
#   bash install.sh    # Standard (user) install
#   bash install.sh 1  # Developer install

platform=$(uname)
if [ "$platform" == 'Linux' ]
then
    profile=~/.bashrc
    miniforge_url=https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
elif [ "$platform" == 'FreeBSD' ] || [ "$platform" == 'Darwin' ]
then
    profile=~/.bash_profile
    miniforge_url=https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-x86_64.sh
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
    echo 'No mamba detected'
    # Is conda installed?
    if ! conda --version
    then
        echo 'No conda detected, either — installing miniforge...'

        # Try to download shell script using curl
        if ! curl -L $miniforge_url -o miniforge.sh
        then
            echo 'Failed to download miniforge installer shell script. Exiting.'
            exit 1
        fi

        # Try to install miniforge
        echo 'Install directory: ~/miniforge'
        if ! bash miniforge.sh -f -b -p ~/miniforge
        then
            echo 'Failed to run miniforge installer shell script. Exiting.'
            exit 1
        fi

        # Set up "activate" command, see:
        # https://docs.anaconda.com/anaconda/install/silent-mode/#linux-macos
        eval "$(~/miniforge/bin/conda shell.bash hook)"
        mamba init

        # Don't automatically activate the environment
        mamba config --set auto_activate_base false

        echo 'mamba installed'
        CONDA_CMD=mamba
    else
        echo 'conda detected'
        CONDA_CMD=conda
    fi
else
    echo 'mamba detected'
    CONDA_CMD=mamba
fi

# This is needed to ensure that environment activation works
source $profile

# Try to activate the base environment
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
    'notebook'
    'obspy'
    'pyqt'
    'rioxarray'
)

# Additional developer packages:
DEVELOPER_PACKAGES=(
    'black=23.10.0'
    'ipython'
    'isort=5.12.0'
    'nbdime'
    'pytest-cov'
    'pytest-mpl'
    'recommonmark'
    'sphinx'
    'sphinx_rtd_theme<3'
    'sphinxcontrib-apidoc'
    'spyder'
    'versioneer'
)

# If user supplied the developer flag, add developer packages to package list
if [ "$1" == 1 ]
then
    PACKAGE_LIST=( "${PACKAGE_LIST[@]}" "${DEVELOPER_PACKAGES[@]}" )
    echo 'Installing developer packages:'
    echo "${DEVELOPER_PACKAGES[@]}"
fi

# Try to create the environment
echo "Creating the $ENV_NAME environment"
if ! $CONDA_CMD create --yes --name $ENV_NAME --channel conda-forge "${PACKAGE_LIST[@]}"
then
    echo 'Failed to create environment. Resolve any conflicts, then try again.'
    exit 1
fi

# Try to activate the new environment
echo "Activating the $ENV_NAME environment"
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

# Try to install this package
echo
echo "Installing $PACKAGE_NAME"
if ! pip install --editable .
then
    echo 'Failed to pip install this package. Exiting.'
    exit 1
fi

# Do additional things if this is a developer install
if [ "$1" == 1 ]
then
    # Give nbdime setup instructions, see:
    # https://nbdime.readthedocs.io/en/latest/#git-integration-quickstart
    echo 'To set up git integration for nbdime, run:'
    echo 'nbdime config-git --enable --global'

    # Install packages only available via pip
    if ! pip install sphinx-markdown-builder
    then
      echo 'Failed to install development pip packages. Exiting.'
      exit 1
    fi
fi
