#!/usr/bin/env bash

# Make sure we fail in case of errors
set -e

# Copy sources (we do not have write permission on the mounted $TRAVIS_BUILD_DIR),
# so let's make a copy of the source code
cd ~
rm -rf my_work_dir
mkdir my_work_dir
# Copy also dot files (.*)
shopt -s dotglob
cp -R ${TRAVIS_BUILD_DIR}/* my_work_dir/

cd my_work_dir

#### borrowed from conda

# Environment
#NUMPYVER=1.15
#MATPLOTLIBVER=2
READLINE_VERSION="6.2"
libgfortranver="3.0"
UPDATE_CONDA=true

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
#python ci/set_minor_version.py --patch $TRAVIS_BUILD_NUMBER --version_file astromodels/version.py

export PKG_VERSION=$(python -c "import versioneer;print(versioneer.get_version())")

echo "HOME= ${HOME}"
echo "Building ${PKG_VERSION} ..."
echo "Python version: ${TRAVIS_PYTHON_VERSION}"
echo "Testing with XSPEC: ${TEST_WITH_XSPEC} ..."

if ${TEST_WITH_XSPEC}; then
    XSPECVER="6.22.1"
    xspec_channel=xspecmodels
    conda config --add channels ${xspec_channel}
    export XSPEC="xspec-modelsonly=${XSPECVER} ${xorg}"
fi

LIBGFORTRAN="libgfortran=${libgfortranver}"

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

echo "dependencies: ${MATPLOTLIB} ${NUMPY} ${XSPEC} ${READLINE} ${LIBGFORTRAN}"


# newer conda is failing hard
#conda install --yes conda=4.5.12

# Answer yes to all questions (non-interactive)
conda config --set always_yes true

# We will upload explicitly at the end, if successful
conda config --set anaconda_upload no

# Create test environment
echo "Create test environment..."
conda create --name test_env -c conda-forge python=$TRAVIS_PYTHON_VERSION pytest codecov pytest-cov git ${MATPLOTLIB} ${NUMPY} ${XSPEC} astropy ${compilers}\
  ${LIBGFORTRAN} scipy pytables krb5=1.14.6 ${READLINE} future

# Make sure conda-forge is the first channel
conda config --add channels defaults

conda config --add channels conda-forge/label/cf201901

conda config --add channels conda-forge

# Activate test environment
echo "Activate test environment..."

source activate test_env

conda config --show channels

conda config --show-sources

# Build package
echo "Build package..."
if $TEST_WITH_XSPEC ; then
    echo "Building WITH xspec"
    if [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
        conda build --python=$TRAVIS_PYTHON_VERSION conda-dist/recipe
        conda index /home/travis/miniconda/conda-bld
    else
    	# there is some strange error about the prefix length
        conda build --no-build-id --python=$TRAVIS_PYTHON_VERSION conda-dist/recipe
        conda index /Users/travis/miniconda/conda-bld
    fi
else
    if [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
	    conda build --python=$TRAVIS_PYTHON_VERSION conda-dist/no_xspec_recipe
	    conda index /home/travis/miniconda/conda-bld
    else
    	# there is some strange error about the prefix length
	    conda build --no-build-id  --python=$TRAVIS_PYTHON_VERSION conda-dist/no_xspec_recipe
	    conda index /Users/travis/miniconda/conda-bld
    fi
fi

echo "======>  installing..."
conda install --use-local -c conda-forge astromodels

if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
    if ${TEST_WITH_XSPEC}; then
        ls /Users/travis/miniconda/envs/test_env/lib/libCCfits*
        ls /Users/travis/miniconda/envs/test_env/lib/libcfitsio*
        ls /Users/travis/miniconda/envs/test_env/lib/libwcs*
        #conda install -c conda-forge/label/cf201901 ccfits=2.5
        ln -s /Users/travis/miniconda/envs/test_env/lib/libCCfits.0.dylib /Users/travis/miniconda/envs/test_env/lib/libCCfits.2.5.dylib
        ls /Users/travis/miniconda/envs/test_env/lib/libCCfits*
    fi
fi

# Run tests
cd astromodels/tests
if ${TEST_WITH_XSPEC}; then
    echo "======>  importing XSPEC..."
    #pytest -s --disable-warnings test_load_xspec_models.py
    python -c "import astromodels.xspec"
fi
python -m pytest -vv --disable-warnings --cov=astromodels # -k "not slow"

# Codecov needs to run in the main git repo

# Upload coverage measurements if we are on Linux
if [[ "$TRAVIS_OS_NAME" == "linux" ]]; then

    echo "********************************** COVERAGE ******************************"
    codecov -t 493c9a2d-42fc-40d6-8e65-24e681efaa1e

fi

# We do not want to upload if we do not test with xspec
if $TEST_WITH_XSPEC ; then
    #echo "======>  importing XSPEC..."
    #python -c "import astromodels.xspec"
    
    # If we are on the master branch upload to the channel
    if [[ "${TRAVIS_EVENT_TYPE}" == "pull_request" ]]; then
        echo "This is a pull request, not uploading to Conda channel"
    
    elif [[ "${TRAVIS_EVENT_TYPE}" == "api" ]]; then
        echo "This build was triggered via API"
    
    elif [[ "${TRAVIS_EVENT_TYPE}" == "push" ]]; then
        echo "This is a push to branch ${TRAVIS_BRANCH}"
        
        echo "${TRAVIS_TAG}"
        if [ -n "${TRAVIS_TAG}" ]; then
            echo "This is the tag ${TRAVIS_TAG}"
            conda install -c conda-forge anaconda-client
            echo "Uploading ${CONDA_BUILD_PATH}"
            
            if [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
                anaconda -v --show-traceback -t $CONDA_UPLOAD_TOKEN upload -u threeml /home/travis/miniconda/conda-bld/linux-64/*.tar.bz2 --force
            else
                anaconda -v --show-traceback -t $CONDA_UPLOAD_TOKEN upload -u threeml /Users/travis/miniconda/conda-bld/*/*.tar.bz2 --force
            fi
        fi
    fi
else
    echo "We didn't test with xspec, not uploading"
    
fi
