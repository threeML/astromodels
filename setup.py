#!/usr/bin/env python

import ctypes.util
import glob
import os
import re
import sys
import subprocess

from packaging import version as packaging_version
from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext as _build_ext

import versioneer

# This is needed to use numpy in this module, and should work whether or not numpy is
# already installed. If it's not, it will trigger an installation

_default_xspec_version = "12.15.0"  # default when installing xspec according following
# https://heasarc.gsfc.nasa.gov/docs/software/conda.html


class My_build_ext(_build_ext):

    def finalize_options(self):

        _build_ext.finalize_options(self)

        import numpy

        self.include_dirs.append(numpy.get_include())
        self.include_dirs.append("astromodels/xspec/include")
        conda_prefix = os.environ.get("CONDA_PREFIX")
        if conda_prefix is not None:
            conda_include_path = os.path.join(conda_prefix, "include")
            self.include_dirs.append(conda_include_path)


def sanitize_lib_name(library_path):
    """Get a fully-qualified library name, like /usr/lib/libgfortran.so.3.0,
    and returns the lib name needed to be passed to the linker in the -l option
    (for example gfortran)

    :param library_path:
    :return:
    """

    lib_name = os.path.basename(library_path)

    # Some regexp magic needed to extract in a system-independent (mac/linux)
    # way the library name

    tokens = re.findall(r"lib(.+)(\.so|\.dylib|\.a)(.+)?", lib_name)

    if not tokens:
        msg = f"Attempting to find {lib_name} in directory {library_path}"
        msg += " but there are no libraries in this directory"

        raise RuntimeError(msg)

    return tokens[0][0]


def find_library(library_root, additional_places=None):
    """Returns the name of the library without extension.

    :param library_root: root of the library to search, for example
        "cfitsio_" will match libcfitsio_1.2.3.4.so
    :return: the name of the library found (NOTE: this is *not* the
        path), and a directory path if the library is not in the system
        paths (and None otherwise). The name of libcfitsio_1.2.3.4.so
        will be cfitsio_1.2.3.4, in other words, it will be what is
        needed to be passed to the linker during a c/c++ compilation, in
        the -l option
    """

    # find_library searches for all system paths in a system independent way (but NOT
    # those defined in LD_LIBRARY_PATH or DYLD_LIBRARY_PATH)

    first_guess = ctypes.util.find_library(library_root)

    if first_guess is not None and library_root == "gfortran":

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
        # (if they sanitize_lib_name(first_guess), are not defined, possible_locations
        # will become [""] which will be handled by the next loop)

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
                # This can happen if there are more than one :, or if nor
                # LD_LIBRARY_PATH nor DYLD_LIBRARY_PATH are defined (because of
                # the default use above for os.environ.get)

                continue

            results = glob.glob(os.path.join(search_path, f"lib{library_root}*"))

            if len(results) >= 1:

                # Results contain things like libXS.so, libXSPlot.so, libXSpippo.so
                # If we are looking for libXS.so, we need to make sure that we get the
                # right one!

                for result in results:

                    if (
                        re.match(
                            f"lib{library_root}" + r"[\-_\.]([0-9])*\d*(\.[0-9]\d*)*",
                            os.path.basename(result),
                        )
                        is None
                    ):

                        continue

                    else:

                        # FOUND IT

                        # This is the full path of the library, like
                        # /usr/lib/libcfitsio_1.2.3.4

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

            # Sanitize the library name to get from the fully-qualified path to just
            # the library name (/usr/lib/libgfortran.so.3.0 becomes gfortran)

            return sanitize_lib_name(library_name), library_dir


def get_xspec_conda_version():
    """Get the version string from conda"""
    try:
        lines = subprocess.check_output(
            ['conda', 'list', '-f', 'xspec']
        ).decode().split('\n')
    except subprocess.CalledProcessError:
        lines = subprocess.check_output(
            ['conda', 'list', '-f', 'xspec']
        ).split('\n')
    for l in lines:
        if not l:
            continue
        if l[0] == '#':
            continue
        tokens = l.split()
        return tokens[1]
    return None


def setup_xspec():

    skip_xspec = os.environ.get("SKIP_XSPEC")
    headas_root = os.environ.get("HEADAS")
    conda_prefix = os.environ.get("CONDA_PREFIX")

    if skip_xspec is not None:

        print(
            "The SKIP_XSPEC env variable was set. Xspec support will not be installed."
        )
        return None

    if headas_root is None:

        # See, maybe we are running in Conda

        if conda_prefix is None:

            # Maybe this is a Conda build

            conda_prefix = os.environ.get("PREFIX")

        if conda_prefix is not None:

            # Yes, this is Conda
            # Let's see if the package xspec-modelsonly has been installed by checking
            # whether one of the Xspec libraries exists within conda
            conda_lib_path = os.path.join(conda_prefix, "lib")
            this_lib, this_lib_path = find_library(
                "XSFunctions", additional_places=[conda_lib_path]
            )

            if this_lib is None:

                # No, there is no library in Conda
                msg = "No xspec-modelsonly package has been installed in Conda. Xspec"
                msg += " support will not be installed"
                print(msg)

                print("Was looking into %s" % conda_lib_path)

                return None

            else:
                msg = ("WARN: The xspec-modelsonly package has been installed"
                       " in Conda, but it's no longer supported."
                       " Xspec support will not be installed")
                print(msg)

                return None

                # Set up the HEADAS variable so that the following will find the
                # libraries
                # headas_root = conda_prefix

        else:

            print("No HEADAS env. variable set. Xspec support will not be installed ")

            return None

    print("HEADAS env. variable detected. Will compile the Xspec extension.")
    print("NOTICE: If you have issues, manually set the environment variable "
          "XSPEC_INC_PATH to the location of the XSPEC headers")
    msg = "If you are still having issues, unset HEADAS before installing and"
    msg += "contact the support team"
    print(msg)

    xspec_version = get_xspec_conda_version()

    if xspec_version is not None:

        print("Found XSPEC version %s in Conda" % xspec_version)

    else:

        print("No XSPEC installation found in Conda")
        print('Xspec was likely compiled from source.')

        xspec_version = os.environ.get("ASTRO_XSPEC_VERSION")

        if xspec_version is None:
            print("WARN: You have not specified an XSPEC version with the ")
            print("WARN: environment variable ASTRO_XSPEC_VERSION")
            print(f"WARN: we will assume you have {_default_xspec_version}")
            print("If you are using a different version of XSPEC, please set"
                  " the environment variable ASTRO_XSPEC_VERSION to the "
                  "version of XSPEC you are using")

            xspec_version = _default_xspec_version

        else:

            print(f"You have specified the XSPEC version {xspec_version}")

    xspec_version = packaging_version.Version(xspec_version)

    if xspec_version < packaging_version.Version("12.12.0"):
        msg = "WARN: XSPEC version is less than 12.12.0, which is the minimal"
        msg += " supported version for astromodels"
        print(msg)
        return None
    elif xspec_version > packaging_version.Version("12.15.0"):
        msg = "WARN: XSPEC version is greater than 12.15.0, which is the"
        msg += " maximal supportedversion for astromodels"
        print(msg)
        return None

    macros = []
    # I am not sure what the naming of the XSPEC components are,
    # but let's stick with major, minor, and patch.
    for major, minor, patch in [
        (12, 12, 0),
        (12, 12, 1),
        (12, 13, 0),
        (12, 13, 1),
        (12, 14, 0),
        (12, 14, 1),
        (12, 15, 0),
    ]:

        version = "{}.{}.{}".format(major, minor, patch)

        macro = "XSPEC_{}_{}_{}".format(major, minor, patch)

        if xspec_version >= packaging_version.Version(version):
            macros += [(macro, None)]

    print(macros)

    # Make sure these libraries exist and are linkable right now
    # (they need to be in LD_LIBRARY_PATH or DYLD_LIBRARY_PATH or in one of the system
    # paths)

    libraries_root = [
        "XSFunctions",
        "XSModel",
        "XSUtil",
        "XS",
        "cfitsio",
        "CCfits",
        "wcs",
        "gfortran",
    ]

    libraries = []
    library_dirs = []

    for lib_root in libraries_root:

        this_library, this_library_path = find_library(
            lib_root, additional_places=[os.path.join(headas_root, "lib")]
        )

        if this_library is None:

            raise IOError(
                "Could not find library %s. Impossible to compile Xspec" % lib_root
            )

        else:

            print("Found library %s in %s" % (this_library, this_library_path))

            libraries.append(this_library)

            if this_library_path is not None:
                # This library is not in one of the system path library, we need to add
                # it to the -L flag during linking. Let's put it in the library_dirs
                # list which will be used in the Extension class

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

        conda_include_path = os.path.join(conda_prefix, "include")
        header_paths.append(conda_include_path)

    # check if there are user set the location of the xspec headers:

    xspec_headers_path = os.environ.get("XSPEC_INC_PATH")

    if xspec_headers_path is not None:

        print("You have set XSPEC_INC_PATH=%s" % xspec_headers_path)

        header_paths.append(xspec_headers_path)

    # Remove duplicates from library_dirs and header_paths

    library_dirs = list(set(library_dirs))
    header_paths = list(set(header_paths))

    print("header paths:")
    for h in header_paths:

        print(f"{h}")

    # Configure the variables to build the external module with the C/C++ wrapper

    ext_modules_configuration = [
        Extension(
            "astromodels.xspec._xspec",
            [
                "astromodels/xspec/src/_xspec.cc",
            ],
            include_dirs=header_paths,
            libraries=libraries,
            library_dirs=library_dirs,
            runtime_library_dirs=library_dirs,
            extra_compile_args=[],
            define_macros=macros,
        ),
    ]

    return ext_modules_configuration


# Normal packages

packages = [
    "astromodels",
    "astromodels/core",
    "astromodels/functions",
    "astromodels/functions/functions_1D",
    "astromodels/functions/dark_matter",
    "astromodels/sources",
    "astromodels/utils",
    "astromodels/xspec",
    "astromodels/tests",
]

# Check whether we can compile Xspec support
ext_modules_configuration = setup_xspec()

# Add the node_ctype module


setup(
    cmdclass=versioneer.get_cmdclass({"build_ext": My_build_ext}),
    version=versioneer.get_version(),
    ext_modules=ext_modules_configuration,
)
