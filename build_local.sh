#!/usr/bin/env bash
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

echo "Running on ${TRAVIS_OS_NAME}"

TEST_WITH_XSPEC=true
TRAVIS_PYTHON_VERSION=2.7
TRAVIS_BUILD_NUMBER=3
ENVNAME=astromodels_test_$TRAVIS_PYTHON_VERSION

# Make sure we fail in case of errors
set -e

# Environment
libgfortranver="3.0"
#XSPECVER="12.10.1b"
#XSPECVER="12.9.1u"
XSPECVER="6.22.1"

#NUMPYVER=1.15
#MATPLOTLIBVER=2
UPDATE_CONDA=true

#xspec_channel=xspec/channel/dev
#xspec_channel=cxc/label/dev

if [[ ${TRAVIS_OS_NAME} == linux ]];
then
    miniconda_os=Linux
    compilers="gcc_linux-64 gxx_linux-64 gfortran_linux-64"
else  # osx
    miniconda_os=MacOSX
    compilers="clang_osx-64 clangxx_osx-64 gfortran_osx-64"

    # On macOS we also need the conda libx11 libraries used to build xspec
    # We also need to pin down ncurses, for now only on macos.
    xorg="xorg-libx11 ncurses=5"
fi




# Get the version in the __version__ environment variable
python ci/set_minor_version.py --patch $TRAVIS_BUILD_NUMBER --version_file astromodels/version.py

export PKG_VERSION=$(cd astromodels && python -c "import version;print(version.__version__)")

echo "Building ${PKG_VERSION} ..."
echo "Python version: ${TRAVIS_PYTHON_VERSION}"
echo "Testing with XSPEC: ${TEST_WITH_XSPEC} ..."

if $UPDATE_CONDA ; then
    # Update conda
    echo "Update conda..."
    conda update --yes -q conda conda-build
fi

#conda config --add channels ${xspec_channel}

if [[ ${TRAVIS_OS_NAME} == osx ]];
then
    conda config --add channels conda-forge
fi

# Figure out requested dependencies
if [ -n "${MATPLOTLIBVER}" ]; then MATPLOTLIB="matplotlib=${MATPLOTLIBVER}"; fi
if [ -n "${NUMPYVER}" ]; then NUMPY="numpy=${NUMPYVER}"; fi
if [ -n "${XSPECVER}" ];
 then export XSPEC="xspec-modelsonly=${XSPECVER} ${xorg}";
fi

echo "dependencies: ${MATPLOTLIB} ${NUMPY}  ${XSPEC}"

# Answer yes to all questions (non-interactive)
conda config --set always_yes true

# We will upload explicitly at the end, if successful
conda config --set anaconda_upload no

# Create test environment
echo "Create test environment..."
#conda create --yes --name $ENVNAME -c conda-forge/label/cf201901 -c fermi python=$TRAVIS_PYTHON_VERSION pytest codecov pytest-cov git ${MATPLOTLIB} ${NUMPY} ${XSPEC} ${compilers}\
#  libgfortran=${libgfortranver} fermitools fermipy

conda create --yes --name $ENVNAME -c conda-forge -c threeml python=$TRAVIS_PYTHON_VERSION pytest codecov pytest-cov git ${MATPLOTLIB} ${NUMPY} ${XSPEC} ${compilers}\
  libgfortran=${libgfortranver}


# Make sure conda-forge is the first channel
conda config --add channels conda-forge

conda config --add channels defaults

# Activate test environment
echo "Activate test environment..."

#source $HOME/work/fermi/miniconda3/etc/profile.d/conda.sh
#conda activate $ENVNAME
source activate $ENVNAME

# Build package
echo "Build package..."
if $TEST_WITH_XSPEC ; then
    echo "Building WITH xspec"
    if [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
        conda build --python=$TRAVIS_PYTHON_VERSION conda-dist/recipe
        #conda index $HOME/work/fermi/miniconda3/conda-bld
        conda index $HOME/miniconda/conda-bld
    else
    	# there is some strange error about the prefix length
        conda build --no-build-id --python=$TRAVIS_PYTHON_VERSION conda-dist/recipe
        conda index $HOME/miniconda/conda-bld
    fi
	echo "======> installing..."
    conda install --use-local -c conda-forge astromodels
    # xspec-modelsonly
else

    if [[ "$TRAVIS_OS_NAME" == "linux" ]]; then

	    conda build --python=$TRAVIS_PYTHON_VERSION conda-dist/no_xspec_recipe
        conda index $HOME/miniconda/conda-bld

    else

	# there is some strange error about the prefix length

	    conda build --no-build-id  --python=$TRAVIS_PYTHON_VERSION conda-dist/no_xspec_recipe
        conda index $HOME/miniconda/conda-bld

    fi

	echo "======>  installing..."
    conda install --use-local -c conda-forge -c threeml astromodels
fi


# Run tests
cd astromodels/tests
python -m pytest -vv --cov=astromodels # -k "not slow"

# Codecov needs to run in the main git repo

# Upload coverage measurements if we are on Linux
if [[ "$TRAVIS_OS_NAME" == "linux" ]]; then

    echo "********************************** COVERAGE ******************************"
    codecov -t 493c9a2d-42fc-40d6-8e65-24e681efaa1e

fi
