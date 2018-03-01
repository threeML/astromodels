#!/usr/bin/env bash

set -e

conda update --yes -q conda
conda config --set always_yes true
conda config --set anaconda_upload no

# Create test environment
conda create --name test_env python=$TRAVIS_PYTHON_VERSION pytest codecov

# Run test
source activate test_env
# Build package
cd /travis_build_dir
conda build -c conda-forge -c threeml --python=$TRAVIS_PYTHON_VERSION conda-dist/recipe
# Install it
conda install --use-local astromodels
# Run tests
cd ~
python -m pytest -vv --cov=astromodels --pyargs astromodels

# Upload coverage measurements
codecov -t 493c9a2d-42fc-40d6-8e65-24e681efaa1e

# If we are on the master branch upload to the channel
if [[ "$TRAVIS_BRANCH" == "linux" ]]; then
        anaconda -t $CONDA_UPLOAD_TOKEN upload -u threeml $CONDA_BLD_PATH/$TRAVIS_OS_NAME-64/*.tar.bz2 --force
else
        echo "On a branch, not uploading to Conda channel"
fi
