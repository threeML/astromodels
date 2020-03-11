#!/usr/bin/env bash
# Make sure we fail in case of errors
set -e

TRAVIS_OS_NAME="unknown"

if [[ "$OSTYPE" == "linux-gnu" ]]; then

        # Linux

        TRAVIS_OS_NAME="linux"


elif [[ "$OSTYPE" == darwin* ]]; then

        # Mac OSX

        TRAVIS_OS_NAME="osx"


elif [[ "$OSTYPE" == "cygwin" ]]; then

        # POSIX compatibility layer and Linux environment emulation for Windows

        TRAVIS_OS_NAME="linux"

else

        # Unknown.

        echo "Could not guess your OS. Exiting."

        exit 1

fi

echo " ===> Running on ${TRAVIS_OS_NAME}"

TEST_WITH_XSPEC=true
USE_LOCAL=false
TRAVIS_PYTHON_VERSION=3.5
TRAVIS_BUILD_NUMBER=6
ENVNAME=astromodels_test_$TRAVIS_PYTHON_VERSION

# Environment
libgfortranver="3.0"

#NUMPYVER=1.15
#MATPLOTLIBVER=2
READLINE_VERSION="6.2"
UPDATE_CONDA=false

if [[ ${TRAVIS_OS_NAME} == linux ]];
then
    miniconda_os=Linux
    compilers="gcc_linux-64 gxx_linux-64 gfortran_linux-64"
else  # osx
    miniconda_os=MacOSX
    compilers="clang_osx-64 clangxx_osx-64 gfortran_osx-64"

    # On macOS we also need the conda libx11 libraries used to build xspec
    # We also need to pin down ncurses, for now only on macos.
    xorg="xorg-libx11" # ncurses=5
fi

# Get the version in the __version__ environment variable
# python ci/set_minor_version.py --patch $TRAVIS_BUILD_NUMBER --version_file astromodels/version.py

# export PKG_VERSION=$(cd astromodels && python -c "import version;print(version.__version__)")


export PKG_VERSION=$(python -c "import versioneer;print(versioneer.get_version())")

echo "Building ${PKG_VERSION} ..."
echo "Python version: ${TRAVIS_PYTHON_VERSION}"
echo "Testing with XSPEC: ${TEST_WITH_XSPEC} ..."
echo "Use local is: ${USE_LOCAL}"

if ${TEST_WITH_XSPEC}; then
    XSPECVER="6.22.1"
    export XSPEC="xspec-modelsonly=${XSPECVER} ${xorg}"
    xspec_channel=xspecmodels
    
    if ${USE_LOCAL}; then
        #conda config --remove channels ${xspec_channel}
        use_local="--use-local"
    else
        conda config --add channels ${xspec_channel}
    fi
fi

if $UPDATE_CONDA ; then
    # Update conda
    echo "Update conda..."
    conda update --yes -q conda conda-build
fi

if [[ ${TRAVIS_PYTHON_VERSION} == 2.7 ]];
then
    READLINE="readline=${READLINE_VERSION}"
fi

# Figure out requested dependencies
if [ -n "${MATPLOTLIBVER}" ]; then MATPLOTLIB="matplotlib=${MATPLOTLIBVER}"; fi
if [ -n "${NUMPYVER}" ]; then NUMPY="numpy=${NUMPYVER}"; fi

echo "dependencies: ${MATPLOTLIB} ${NUMPY}  ${XSPEC}"

# Answer yes to all questions (non-interactive)
conda config --set always_yes true

# We will upload explicitly at the end, if successful
conda config --set anaconda_upload no

# Create test environment
echo "Create test environment..."

conda create --yes --name $ENVNAME -c conda-forge ${use_local} python=$TRAVIS_PYTHON_VERSION pytest codecov pytest-cov git ${MATPLOTLIB} ${NUMPY} ${XSPEC} astropy ${compilers}\
  libgfortran=${libgfortranver} scipy pytables krb5=1.14.6 ${READLINE} future


# Make sure conda-forge is the first channel
conda config --add channels defaults

conda config --add channels conda-forge/label/cf201901

conda config --add channels conda-forge

# Activate test environment
echo "Activate test environment..."

source $CONDA_PREFIX/etc/profile.d/conda.sh
#source /home/ndilalla/work/fermi/miniconda3/etc/profile.d/conda.sh
conda activate $ENVNAME

# Build package
echo "Build package..."
conda config --show channels

if $TEST_WITH_XSPEC ; then
    echo " ====> Building WITH xspec"
    if [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
        conda build --python=$TRAVIS_PYTHON_VERSION conda-dist/recipe
    else
    	# there is some strange error about the prefix length
        conda build --no-build-id --python=$TRAVIS_PYTHON_VERSION conda-dist/recipe
        #conda install -c conda-forge/label/cf201901 ccfits=2.5
    fi
else
    echo " ====> Building WITHOUT xspec"
    if [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
	    conda build --python=$TRAVIS_PYTHON_VERSION conda-dist/no_xspec_recipe
    else
    	# there is some strange error about the prefix length
	    conda build --no-build-id  --python=$TRAVIS_PYTHON_VERSION conda-dist/no_xspec_recipe
    fi
fi

echo "======>  installing..."
conda install --use-local -c conda-forge astromodels

echo "======>  Run tests..."

cd astromodels/tests
python -m pytest -vv --cov=astromodels # -k "not slow"

# Codecov needs to run in the main git repo

# Upload coverage measurements if we are on Linux
if [[ "$TRAVIS_OS_NAME" == "linux" ]]; then

    echo "********************************** COVERAGE ******************************"
    codecov -t 493c9a2d-42fc-40d6-8e65-24e681efaa1e

fi
