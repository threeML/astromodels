#!/usr/bin/env bash

# This ensure that the script will exit if any command fails
set -e

# I am running as root, first let's create a new user and become that
# The user_id env. variable must be specified on the docker command line using -e user_id=`id -u`
adduser --system --home /home/user --shell /bin/bash --uid $user_id user --disabled-password
exec sudo -i -u user /bin/bash << EOF

###################################################################################
# Beginning of script run as user "user"
###################################################################################

set -e
cd /home/user

# Setup environment
echo "##########################################################"
echo " Setting up HEASOFT environment"
echo "##########################################################"

export HEADAS=/heasoft/build/x86_64-unknown-linux-gnu-libc2.23-0
source /heasoft/build/x86_64-unknown-linux-gnu-libc2.23-0/headas-init.sh

# Print XSPEC version
echo exit | xspec

echo "##########################################################"
echo " Creating python virtual environment"
echo "##########################################################"
virtualenv astromodels_env
source astromodels_env/bin/activate

echo "##########################################################"
echo " Installing numpy, pytest, pytest-cov and coveralls"
echo "##########################################################"

pip install numpy pytest pytest-cov coveralls codecov

echo "##########################################################"
echo " Installing astromodels"
echo "##########################################################"

# This assumes that the $TRAVIS_BUILD_DIR directory has been mounted to
# /astromodels using -v $TRAVIS_BUILD_DIR:/travis_build_dir
cd /travis_build_dir
pip install .

echo "##########################################################"
echo " Executing tests and coveralls"
echo "##########################################################"

# Execute tests
# (need to move away from root directory, otherwise xspec import will fail)
cd astromodels
python -m pytest -vv --cov=astromodels

echo "##########################################################"
echo " Executing codecov"
echo "##########################################################"

# Execute the coverage analysis
codecov -t 493c9a2d-42fc-40d6-8e65-24e681efaa1e

###################################################################################
# end of script run as user "user"
###################################################################################


EOF