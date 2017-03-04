#!/usr/bin/env python

import ctypes.util
import glob
import sys

import os
import re
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext as _build_ext


# This is needed to use numpy in this module, and should work whether or not numpy is
# already installed. If it's not, it will trigger an installation

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

    return tokens[0][0]


def find_library(library_root):
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


# Get the version number
execfile('astromodels/version.py')


def setup_xspec():
    headas_root = os.environ.get("HEADAS")

    if headas_root is None:

        print("No HEADAS env. variable set. Xspec support will not be installed ")

        return None

    else:

        print("\n Xspec is detected. Will compile the Xspec extension.\n")

    # Make sure these libraries exist and are linkable right now
    # (they need to be in LD_LIBRARY_PATH or DYLD_LIBRARY_PATH or in one of the system paths)

    libraries_root = ['XSFunctions', 'XSModel', 'XSUtil', 'XS', 'cfitsio', 'CCfits', 'wcs', 'gfortran']

    libraries = []
    library_dirs = []

    for lib_root in libraries_root:

        this_library, this_library_path = find_library(lib_root)

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

    # Remove duplicates from library_dirs

    library_dirs = list(set(library_dirs))

    # Configure the variables to build the external module with the C/C++ wrapper
    ext_modules_configuration = [

        Extension("astromodels.xspec._xspec",

                  ["astromodels/xspec/src/_xspec.cc", ],

                  libraries=libraries,
                  library_dirs=library_dirs,
                  extra_compile_args=[])]

    return ext_modules_configuration


# Normal packages

packages = ['astromodels',
            'astromodels/core',
            'astromodels/functions',
            'astromodels/functions/dark_matter',
            'astromodels/sources',
            'astromodels/utils',
            'astromodels/xspec',
            'astromodels/tests'
            ]

# Check whether we can compile Xspec support
ext_modules_configuration = setup_xspec()

# Add the node_ctype module

# This defines the external module
node_ctype_ext = Extension('astromodels.core.node_ctype',
                           sources = ['astromodels/core/node_ctype/node_ctype.cxx'],
                           extra_compile_args=[]) # '-UNDEBUG' for debugging


if ext_modules_configuration is None:

    # No Xspec
    ext_modules_configuration = [node_ctype_ext]

else:

    ext_modules_configuration.append(node_ctype_ext)


setup(
    name="astromodels",

    setup_requires=['numpy'],

    cmdclass={'build_ext': My_build_ext},

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
        'numdifftools',
        'tables',
        'pandas',
        'html2text'],

    ext_modules=ext_modules_configuration,

    package_data={
        'astromodels': ['data/dark_matter/*'],
    },
    include_package_data=True,

)
