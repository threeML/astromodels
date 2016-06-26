#!/usr/bin/env python

import os
import sys
import glob
import numpy as np
import re

from setuptools import setup, Extension

# Get the version number
execfile('astromodels/version.py')


def setup_xspec():

    headas_root = os.environ.get("HEADAS")

    if headas_root is None:

        print("No HEADAS env. variable set. Xspec support will not be installed ")

        return None

    library_dirs = [os.path.join(headas_root, 'lib')]

    # Check that the library directories exist

    if not os.path.exists(library_dirs[0]):

        print("\nERROR: the library directory %s of HEADAS does not exist!" % library_dirs[0])

        sys.exit(-1)

    libraries = ['XSFunctions', 'XSModel', 'XSUtil', 'XS']

    # Check that the libraries exist
    for library in libraries:
        
        # Check for Linux/UNIX
        
        library_name = 'lib%s.so' % library

        if not os.path.exists(os.path.join(library_dirs[0], library_name)):

            # See if it has been compiled static (usually libXS is)

            library_name = 'lib%s.a' % library

            if not os.path.exists(os.path.join(library_dirs[0], library_name)):

                # If on OS X look for dylib
                
                library_name = 'lib%s.dylib' % library

                if not os.path.exists(os.path.join(library_dirs[0], library_name)):


                    print("\nERROR: the library %s does not exist in %s while setting up Xspec!" %
                          (library, library_dirs[0]))

    # Now find versions for required "external" libraries (part of HEASOFT but not Xspec by itself)

    required_libraries = ['cfitsio_', 'CCfits_', 'wcs-']
    

    for library_to_probe in required_libraries:
        
        # Linux/UNIX
        
        search_path = os.path.join(library_dirs[0], 'lib%s*.so' % library_to_probe)
        
        versions = glob.glob(search_path)

        if len(versions) == 0:
            
            search_path = os.path.join(library_dirs[0], 'lib%s*.so' % (library_to_probe[:-1]+'.'))
            
            versions = glob.glob(search_path)

            if len(versions) == 0:
                
                # Mac
                
                search_path = os.path.join(library_dirs[0], 'lib%s*.dylib' % library_to_probe)
                
                versions = glob.glob(search_path)

                if len(versions) == 0:

                    search_path = os.path.join(library_dirs[0], 'lib%s*.dylib' % (library_to_probe[:-1]+'.'))
                    
                    versions = glob.glob(search_path)

                    if len(versions) == 0:

                        print("\nERROR: cannot find version for library %s while setting up Xspec" % (library_to_probe))
                        sys.exit(-1)

        # Up to there versions[0] is a fully-qualified path
        # we need instead just the name of the library, without
        # the lib prefix nor the .so nor .a nor .dylib extension

        name = os.path.basename(versions[0])
        sanitized_name = re.match('lib(.+)\.', name).groups()[0]

        libraries.append(sanitized_name)

    # We also need gfortran, which must be installed on his own. If there is Xspec installed,
    # then there is also gfortran somewhere because it is a dependence
    libraries.append('gfortran')

    # Configure the variables to build the external module with the C/C++ wrapper

    ext_modules_configuration = [

        Extension("astromodels.xspec._xspec",

                  ["astromodels/xspec/src/_xspec.cc", ],

                  libraries=libraries,

                  include_dirs=['astromodels/xspec/include', np.get_include()],

                  library_dirs=library_dirs,
                  extra_compile_args=[])]

    return ext_modules_configuration

# Normal packages

packages = ['astromodels',
            'astromodels/functions',
            'astromodels/sources',
            'astromodels/utils',
            'astromodels/xspec'
            ]

# Check whether we can compile Xspec support
ext_modules_configuration = setup_xspec()

setup(
    name="astromodels",

    packages=packages,

    data_files=[('astromodels/data/functions', glob.glob('astromodels/data/functions/*.yaml'))],

    # The __version__ comes from the exec at the top

    version=__version__,

    description="Astromodels contains models to be used in likelihood or Bayesian analysis in astronomy",

    author='Giacomo Vianello',

    author_email='giacomo.vianello@gmail.com',

    url='https://github.com/giacomov/astromodels',

    download_url='https://github.com/giacomov/astromodels/archive/v0.1',

    keywords=['Likelihood', 'Models', 'fit'],

    classifiers=[],

    install_requires=[
        'numpy >= 1.6',
        'PyYAML',
        'astropy >= 1.0',
        'scipy>=0.13',
        'numdifftools'],

    ext_modules=ext_modules_configuration

)
