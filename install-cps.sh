#!/usr/bin/env bash

# USAGE:
#   bash install-cps.sh

# NOTES:
# Computer Programs in Seismology (CPS) is a collection of programs which must be
# compiled. CPS attempts to compile every program. If an error is encountered, the next
# program is attempted. Thus, the compilation step will always run "successfully",
# regardless of how many programs are compiled. At the end of script, the number of
# compiled executables is printed. It is the responsibility of the user to verify that
# all their needed programs exist. If some programs are missing, the user should check
# the file C.txt in the PROGRAMS.330/ directory — searching for "error" should be enough
# to reveal the problem. In almost all cases, problems arise from missing dependencies
# required for compilation. Use the error messages to understand what needs to be
# installed, and/or consult the CPS documentation at:
# https://rbherrmann.github.io/ComputerProgramsSeismology/

# URL pointing to current CPS release
CPS_URL='https://rbherrmann.github.io/ComputerProgramsSeismology/NP330.Mar-29-2025.tgz'

# Target directory to contain CPS directory, `PROGRAMS.330/`
CPS_LOCATION=~

# Determine the platform for CPS compile setup
platform=$(uname)
if [ "$platform" == 'Linux' ]
then
    setup_flag='LINUX6440'
elif [ "$platform" == 'Darwin' ]
then
    if [ "$(uname -m)" == 'arm64' ]
    then
        setup_flag='OSXM'
    else
        setup_flag='OSXIntel'
    fi
else
    echo 'Unsupported platform. Exiting.'
    exit 1
fi

# Check for existing CPS directory
if [ -d $CPS_LOCATION/PROGRAMS.330 ]
then
    printf "CPS directory already exists at: $CPS_LOCATION/PROGRAMS.330/\nRemove it and try again — exiting.\n"
    exit 1
fi

# Temporary tarball location
cps_tarball=$(mktemp)

# Download the CPS tarball
curl --location $CPS_URL --output $cps_tarball

# Try to unzip the CPS tarball
if ! tar --directory $CPS_LOCATION -xzf $cps_tarball
then
    printf "\nFailed to unzip CPS tarball — bad URL? Check:\n$CPS_URL\nExiting.\n"
    rm $cps_tarball
    exit 1
fi
rm $cps_tarball

# Move into the CPS directory
cd $CPS_LOCATION/PROGRAMS.330

# Set up for the proper platform (need to override $SHELL)
SHELL=$(which bash)
./Setup $setup_flag

# Compile
printf '\nCompiling CPS...\n'
./C > C.txt 2>&1
tail -1 C.txt  # The last line of the log file gives the number of programs installed
printf '\nCPS installation finished. Add the following location to your PATH:\n'
echo $CPS_LOCATION/PROGRAMS.330/bin/
