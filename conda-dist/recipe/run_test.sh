#!/usr/bin/bash

arch=$(uname -m)
os=$(uname -o)

if [[ $arch != 'x86_64' && $os == 'GNU/Linux' ]];
then
	echo "Detected Linux aarch64, skippinx XSPEC tests";
else
	conda install -c https://heasarc.gsfc.nasa.gov/FTP/software/conda/ -c conda-forge xspec==12.15.1 -y;
fi;

pytest -vv astromodels/tests --cov=astromodels --cov-report=term;
