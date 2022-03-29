#!/usr/bin/env python

from __future__ import print_function

import ctypes.util
import glob
import os
import re
import sys
from distutils.version import LooseVersion

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext as _build_ext

import versioneer

# This is needed to use numpy in this module, and should work whether or not numpy is
# already installed. If it's not, it will trigger an installation

_default_xspec_version = "12.10.1"


class My_build_ext(_build_ext):

    def finalize_options(self):

        _build_ext.finalize_options(self)

        # Prevent numpy from thinking it is still in its setup process:

        __builtins__.__NUMPY_SETUP__ = False

        import numpy

        self.include_dirs.append(numpy.get_include())
        self.include_dirs.append('astromodels/xspec/include')


def sanitize_lib_name(library_path):
    """
    Get a fully-qualified library name, like /usr/lib/libgfortran.so.3.0, and returns the lib name needed to be
    passed to the linker in the -l option (for example gfortran)

    :param library_path:
    :return:
    """

    lib_name = os.path.basename(library_path)

    # Some regexp magic needed to extract in a system-independent (mac/linux) way the library name

    tokens = re.findall("lib(.+)(\.so|\.dylib|\.a)(.+)?", lib_name)

    if not tokens:
        raise RuntimeError('Attempting to find %s in directory %s but there are no libraries in this directory'%(lib_name,library_path))


    return tokens[0][0]


def find_library(library_root, additional_places=None):
    """
    Returns the name of the library without extension

    :param library_root: root of the library to search, for example "cfitsio_" will match libcfitsio_1.2.3.4.so
    :return: the name of the library found (NOTE: this is *not* the path), and a directory path if the library is not
    in the system paths (and None otherwise). The name of libcfitsio_1.2.3.4.so will be cfitsio_1.2.3.4, in other words,
    it will be what is needed to be passed to the linker during a c/c++ compilation, in the -l option
    """

    # find_library searches for all system paths in a system independent way (but NOT those defined in
    # LD_LIBRARY_PATH or DYLD_LIBRARY_PATH)

    first_guess = ctypes.util.find_library(library_root)

    if first_guess is not None:

        # Found in one of the system paths

        if sys.platform.lower().find("linux") >= 0:

            # On linux the linker already knows about these paths, so we
            # can return None as path

            return sanitize_lib_name(first_guess), None

        elif sys.platform.lower().find("darwin") >= 0:

            # On Mac we still need to return the path, because the linker sometimes
            # does not look into it

            return sanitize_lib_name(first_guess), os.path.dirname(first_guess)

        else:

            # Windows is not supported

            raise NotImplementedError("Platform %s is not supported" % sys.platform)

    else:

        # could not find it. Let's examine LD_LIBRARY_PATH or DYLD_LIBRARY_PATH
        # (if they sanitize_lib_name(first_guess), are not defined, possible_locations will become [""] which will
        # be handled by the next loop)

        if sys.platform.lower().find("linux") >= 0:

            # Unix / linux

            possible_locations = os.environ.get("LD_LIBRARY_PATH", "").split(":")

        elif sys.platform.lower().find("darwin") >= 0:

            # Mac

            possible_locations = os.environ.get("DYLD_LIBRARY_PATH", "").split(":")

        else:

            raise NotImplementedError("Platform %s is not supported" % sys.platform)

        if additional_places is not None:

            possible_locations.extend(additional_places)

        # Now look into the search paths

        library_name = None
        library_dir = None

        for search_path in possible_locations:

            if search_path == "":
                # This can happen if there are more than one :, or if nor LD_LIBRARY_PATH
                # nor DYLD_LIBRARY_PATH are defined (because of the default use above for os.environ.get)

                continue

            results = glob.glob(os.path.join(search_path, "lib%s*" % library_root))

            if len(results) >= 1:

                # Results contain things like libXS.so, libXSPlot.so, libXSpippo.so
                # If we are looking for libXS.so, we need to make sure that we get the right one!

                for result in results:

                    if re.match("lib%s[\-_\.]" % library_root, os.path.basename(result)) is None:

                        continue

                    else:

                        # FOUND IT

                        # This is the full path of the library, like /usr/lib/libcfitsio_1.2.3.4

                        library_name = result
                        library_dir = search_path

                        break

            else:

                continue

            if library_name is not None:
                break

        if library_name is None:

            return None, None

        else:

            # Sanitize the library name to get from the fully-qualified path to just the library name
            # (/usr/lib/libgfortran.so.3.0 becomes gfortran)

            return sanitize_lib_name(library_name), library_dir




def setup_xspec():

    headas_root = os.environ.get("HEADAS")
    conda_prefix = os.environ.get("CONDA_PREFIX")
    xspec_version = os.environ.get("ASTRO_XSPEC_VERSION")


    # thanks to the sherpa team for this
    
    if xspec_version is None:

        print("WARN: You have not specified and XSPEC version with the ")
        print("WARN: environment variable ASTRO_XSPEC_VERSION")
        print(f"WARN: we will assume you have {_default_xspec_version}")

        xspec_raw_version = _default_xspec_version 

    else:

        print(f"WARN: you have specified you have XSPEC version {xspec_version}")

        xspec_raw_version = xspec_version


    xspec_version = LooseVersion(xspec_raw_version)

    macros = []

    if xspec_version < LooseVersion("12.9.0"):
        print("WARN: XSPEC Version is less than 12.9.0, which is the minimal supported version for astromodels")

        # I am not sure what the naming of the XSPEC components are,
        # but let's stick with major, minor, and patch.
        #
    for major, minor, patch in [(12, 9, 0), (12, 9, 1),
                                (12, 10, 0), (12, 10, 1),
                                (12, 11, 0), (12, 11, 1),
                                (12, 12, 0), (12, 12, 1)]:

        version = '{}.{}.{}'.format(major, minor, patch)

        macro = 'XSPEC_{}_{}_{}'.format(major, minor, patch)

        if xspec_version >= LooseVersion(version):
            macros += [(macro, None)]
                        
    print(macros)
                
    if headas_root is None:

        # See, maybe we are running in Conda
        
        if conda_prefix is None:
            
            # Maybe this is a Conda build
            
            conda_prefix = os.environ.get("PREFIX")

        if conda_prefix is not None:

            # Yes, this is Conda
            # Let's see if the package xspec-modelsonly has been installed by checking whether one of the Xspec
            # libraries exists within conda
            conda_lib_path = os.path.join(conda_prefix, 'lib')
            this_lib, this_lib_path = find_library('XSFunctions', additional_places=[conda_lib_path])

            if this_lib is None:

                # No, there is no library in Conda
                print("No xspec-modelsonly package has been installed in Conda. Xspec support will not be installed")

                print("Was looking into %s" % conda_lib_path)

                return None

            else:

                print("The xspec-modelsonly package has been installed in Conda. Xspec support will be installed")

                # Set up the HEADAS variable so that the following will find the libraries
                headas_root = conda_prefix

        else:

            print("No HEADAS env. variable set. Xspec support will not be installed ")

            return None

    else:

        print("\n Xspec is detected. Will compile the Xspec extension.\n")
        print("\n NOTICE!!!!!\n")
        print("If you have issues, manually set the ENV variable XSPEC_INC_PATH")
        print("To the location of the XSPEC headers\n\n")
        print("If you are still having issues, unset HEADAS before installing and contact the support team")
        


    # Make sure these libraries exist and are linkable right now
    # (they need to be in LD_LIBRARY_PATH or DYLD_LIBRARY_PATH or in one of the system paths)
    
    libraries_root = ['XSFunctions', 'XSModel', 'XSUtil', 'XS', 'cfitsio', 'CCfits', 'wcs', 'gfortran']
            
    libraries = []
    library_dirs = []

    for lib_root in libraries_root:

        this_library, this_library_path = find_library(lib_root, additional_places=[os.path.join(headas_root, 'lib')])

        if this_library is None:

            raise IOError("Could not find library %s. Impossible to compile Xspec" % lib_root)

        else:

            print("Found library %s in %s" % (this_library, this_library_path))

            libraries.append(this_library)

            if this_library_path is not None:
                # This library is not in one of the system path library, we need to add
                # it to the -L flag during linking. Let's put it in the library_dirs list
                # which will be used in the Extension class

                library_dirs.append(this_library_path)


    # try to manually add on the include directory

    header_paths = []
    
    if library_dirs:

        # grab it from the lib assuming that it is one up
        xspec_path, _ = os.path.split(library_dirs[0])
        include_path = os.path.join(xspec_path, "include")

        header_paths.append(include_path)

    # let's be sure to add the conda include directory
       
    if conda_prefix is not None:
        
        conda_include_path = os.path.join(conda_prefix, 'include')
        header_paths.append(conda_include_path)

    # check if there are user set the location of the xspec headers:
   
    xspec_headers_path = os.environ.get("XSPEC_INC_PATH")
    
    if xspec_headers_path is not None:

        print("You have set XSPEC_INC_PATH=%s" % xspec_headers_path)

        header_paths.append(xspec_headers_path)
    
    # Remove duplicates from library_dirs and header_paths

    library_dirs = list(set(library_dirs))  
    header_paths = list(set(header_paths))

    # Configure the variables to build the external module with the C/C++ wrapper


    ext_modules_configuration = [

        Extension("astromodels.xspec._xspec",

                  ["astromodels/xspec/src/_xspec.cc", ],
                  include_dirs=header_paths,
                  libraries=libraries,
                  library_dirs=library_dirs,
                  runtime_library_dirs=library_dirs,
                  extra_compile_args=[], define_macros=macros),
    ]

    return ext_modules_configuration


# Normal packages

packages = ['astromodels',
            'astromodels/core',
            'astromodels/functions',
            'astromodels/functions/functions_1D',
            'astromodels/functions/dark_matter',
            'astromodels/sources',
            'astromodels/utils',
            'astromodels/xspec',
            'astromodels/tests'
            ]

# Check whether we can compile Xspec support
ext_modules_configuration = setup_xspec()

# Add the node_ctype module


setup(

    setup_requires=['numpy'],

    #cmdclass={'build_ext': My_build_ext},
    cmdclass=versioneer.get_cmdclass({'build_ext': My_build_ext}),
    
    packages=packages,

    data_files=[('astromodels/data/functions', glob.glob('astromodels/data/functions/*.yaml')),
                ('astromodels/data/tests',  glob.glob('astromodels/data/tests/*.fits'))

    ],

    # The __version__ comes from the exec at the top

    version=versioneer.get_version(),


    download_url='https://github.com/threeml/astromodels/archive/v0.1',

    keywords=['Likelihood', 'Models', 'fit'],
    

    ext_modules=ext_modules_configuration,

    package_data={
        'astromodels': ['data/dark_matter/*', 'data/xsect/*', 'data/past_1D_values.h5'],
    },


)
