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

_default_xspec_version = "12.15.1"  # default when installing xspec according following
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

    def run(self):
        """Build extensions and fix library paths on macOS"""
        _build_ext.run(self)

        # On macOS, fix library paths for extensions that use extra_objects
        if sys.platform.lower().find("darwin") >= 0:
            for ext in self.extensions:
                if ext.name == "astromodels.xspec._xspec" and ext.extra_objects:
                    self._fix_macos_library_paths(ext)

    def _fix_macos_library_paths(self, ext):
        """Fix library install names on macOS using install_name_tool"""
        import subprocess

        # Get the built extension file
        build_lib = self.get_ext_fullpath(ext.name)
        if not os.path.exists(build_lib):
            return

        # Check what libraries the extension depends on using otool
        try:
            otool_output = subprocess.check_output(
                ["otool", "-L", build_lib], stderr=subprocess.STDOUT
            ).decode()
        except (subprocess.CalledProcessError, FileNotFoundError):
            # otool might not be available, skip fixing
            return

        # Parse otool output to find library references
        # otool -L output format:
        # /path/to/lib.dylib (compatibility version X.Y.Z, current version A.B.C)
        #     /path/to/another/lib.dylib (compatibility version ...)
        library_refs = []
        for line in otool_output.split("\n"):
            line = line.strip()
            if line and not line.startswith(build_lib):
                # Extract the library path/name (everything before the first space or
                # parent)
                if "(" in line:
                    lib_ref = line.split("(")[0].strip()
                else:
                    lib_ref = line.strip()
                if lib_ref:
                    library_refs.append(lib_ref)

        print(f"Found {len(library_refs)} library references in extension")

        # For each library in extra_objects, fix references to point to
        # the actual library file that exists
        for lib_path in ext.extra_objects:
            if not os.path.exists(lib_path):
                continue

            lib_dir = os.path.dirname(lib_path)
            lib_name = os.path.basename(lib_path)

            # Check if this library is referenced in the extension
            # Look for references that might point to a versioned library
            # that doesn't exist (e.g., libwcs.8.dylib when only
            # libwcs.8.3.dylib exists)
            match = re.match(r"lib(.+?)\.(\d+)\.(\d+)\.dylib$", lib_name)
            if match:
                base_name = match.group(1)
                major_version = match.group(2)

                # Check what install name the library file itself has
                # This tells us what reference will be embedded in our extension
                try:
                    otool_d_output = subprocess.check_output(
                        ["otool", "-D", lib_path], stderr=subprocess.STDOUT
                    ).decode()
                    # otool -D output format:
                    # lib_path:
                    #     install_name (or multiple if there are multiple)
                    install_names = []
                    for line in otool_d_output.split("\n"):
                        line = line.strip()
                        if line and not line.endswith(":"):
                            install_names.append(line)
                except (subprocess.CalledProcessError, FileNotFoundError):
                    install_names = []

                # The major version library name we're looking for
                major_version_lib = f"lib{base_name}.{major_version}.dylib"
                new_ref = f"@rpath/{lib_name}"

                # Find all references that match the major version library
                # Try to fix them by changing to the actual library file
                fixed = False
                for lib_ref in library_refs:
                    # Check if this reference matches the major version library
                    # (could be full path, relative path, or just the name)
                    if major_version_lib in lib_ref or lib_ref.endswith(
                        major_version_lib
                    ):
                        # Check if the major version library file doesn't exist
                        old_ref = None
                        if lib_ref.startswith("/"):
                            # Full path reference
                            if not os.path.exists(lib_ref):
                                # The referenced file doesn't exist, change it
                                old_ref = lib_ref
                        else:
                            # Relative or just name reference
                            # Check if it exists in the library directory
                            potential_path = os.path.join(lib_dir, lib_ref)
                            if not os.path.exists(potential_path):
                                old_ref = lib_ref

                        if old_ref:
                            # Change the reference to point to the actual library file
                            # using @rpath so it can find it at runtime
                            try:
                                subprocess.check_call(
                                    [
                                        "install_name_tool",
                                        "-change",
                                        old_ref,
                                        new_ref,
                                        build_lib,
                                    ],
                                    stderr=subprocess.DEVNULL,
                                )
                                print(
                                    f"Fixed library reference: {old_ref} -> {new_ref}"
                                )
                                fixed = True
                            except (subprocess.CalledProcessError, FileNotFoundError):
                                # install_name_tool might not be available or might fail
                                print(
                                    f"Warning: Could not fix library reference "
                                    f"{old_ref} (install_name_tool failed)"
                                )

                # If we didn't fix it yet, try a more aggressive approach:
                # change any reference that contains the major version library name
                if not fixed:
                    for lib_ref in library_refs:
                        if major_version_lib in lib_ref:
                            try:
                                subprocess.check_call(
                                    [
                                        "install_name_tool",
                                        "-change",
                                        lib_ref,
                                        new_ref,
                                        build_lib,
                                    ],
                                    stderr=subprocess.DEVNULL,
                                )
                                print(
                                    f"Fixed library reference (fallback): "
                                    f"{lib_ref} -> {new_ref}"
                                )
                                fixed = True
                                break
                            except (subprocess.CalledProcessError, FileNotFoundError):
                                # Try next reference
                                pass


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

    tokens = re.findall(r"lib(.+)(\.so|\.dylib|\.a|\.la)(.+)?", lib_name)

    if not tokens:
        msg = f"Attempting to find {lib_name} in directory {library_path}"
        msg += " but there are no libraries in this directory"

        raise RuntimeError(msg)

    return tokens[0][0]


def find_library(library_root, additional_places=None):
    """Returns the name of the library without extension.

    :param library_root: root of the library to search, for example
        "cfitsio_" will match libcfitsio_1.2.3.4.so
    :return: a tuple of (library_name, library_dir, full_path) where:
        - library_name: the name to be passed to the linker in the
          -l option
        - library_dir: the directory path (None if in system paths)
        - full_path: the full path to the library file (used when
          unversioned symlink is missing)
    """

    # find_library searches for all system paths in a system independent way (but NOT
    # those defined in LD_LIBRARY_PATH or DYLD_LIBRARY_PATH)

    first_guess = ctypes.util.find_library(library_root)

    if first_guess is not None and library_root == "gfortran":

        # Found in one of the system paths

        if sys.platform.lower().find("linux") >= 0:

            # On linux the linker already knows about these paths, so we
            # can return None as path

            return sanitize_lib_name(first_guess), None, None

        elif sys.platform.lower().find("darwin") >= 0:

            # On Mac we still need to return the path, because the linker sometimes
            # does not look into it

            return sanitize_lib_name(first_guess), os.path.dirname(first_guess), None

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
        library_full_path = None

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

                        library_full_path = result
                        library_name = result
                        library_dir = search_path

                        break

            else:

                continue

            if library_name is not None:
                break

        if library_name is None:

            return None, None, None

        else:

            # Sanitize the library name to get from the fully-qualified path
            # to just the library name (/usr/lib/libgfortran.so.3.0 becomes
            # gfortran)

            sanitized_name = sanitize_lib_name(library_name)

            # Check if the unversioned symlink exists
            # If not, we'll need to pass the full path to the linker
            if sys.platform.lower().find("linux") >= 0:
                extension = ".so"
            elif sys.platform.lower().find("darwin") >= 0:
                extension = ".dylib"
            else:
                extension = ".so"

            unversioned_lib = os.path.join(
                library_dir, f"lib{sanitized_name}{extension}"
            )

            if os.path.exists(unversioned_lib):
                # Unversioned symlink exists, use normal linking
                return sanitized_name, library_dir, None
            else:
                # No unversioned symlink exists
                # On macOS, check for major version symlink (e.g., libwcs.8.dylib)
                if sys.platform.lower().find("darwin") >= 0:
                    # Extract major version from library name if present
                    # e.g., libwcs.8.3.dylib -> try libwcs.8.dylib
                    base_name = os.path.basename(library_full_path)
                    # Match pattern like libwcs.8.3.dylib or libwcs.so.8.3
                    version_match = re.search(r"\.(\d+)\.(\d+)\.dylib$", base_name)
                    if version_match:
                        major_version = version_match.group(1)
                        major_version_lib = os.path.join(
                            library_dir, f"lib{sanitized_name}.{major_version}.dylib"
                        )
                        if os.path.exists(major_version_lib):
                            # Use major version symlink instead of full version
                            return sanitized_name, library_dir, major_version_lib

                # No suitable symlink found, return the full path for direct
                # linking
                return sanitized_name, library_dir, library_full_path


def get_xspec_conda_version():
    """Get the version string from conda"""
    try:
        lines = (
            subprocess.check_output(["conda", "list", "-f", "xspec"])
            .decode()
            .split("\n")
        )
    except subprocess.CalledProcessError:
        lines = subprocess.check_output(["conda", "list", "-f", "xspec"]).split("\n")
    for l in lines:
        if not l:
            continue
        if l[0] == "#":
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
            this_lib, this_lib_path, full_lib_library = find_library(
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
                msg = (
                    "WARN: The xspec-modelsonly package has been installed"
                    " in Conda, but it's no longer supported."
                    " Xspec support will not be installed"
                )
                print(msg)

                return None

                # Set up the HEADAS variable so that the following will find the
                # libraries
                # headas_root = conda_prefix

        else:

            print("No HEADAS env. variable set. Xspec support will not be installed ")

            return None

    print("HEADAS env. variable detected. Will compile the Xspec extension.")
    print(
        "NOTICE: If you have issues, manually set the environment variable "
        "XSPEC_INC_PATH to the location of the XSPEC headers"
    )
    msg = "If you are still having issues, unset HEADAS before installing and"
    msg += "contact the support team"
    print(msg)

    xspec_version = get_xspec_conda_version()

    if xspec_version is not None:

        print("Found XSPEC version %s in Conda" % xspec_version)

    else:

        print("No XSPEC installation found in Conda")
        print("Xspec was likely compiled from source.")

        xspec_version = os.environ.get("ASTRO_XSPEC_VERSION")

        if xspec_version is None:
            print("WARN: You have not specified an XSPEC version with the ")
            print("WARN: environment variable ASTRO_XSPEC_VERSION")
            print(f"WARN: we will assume you have {_default_xspec_version}")
            print(
                "If you are using a different version of XSPEC, please set"
                " the environment variable ASTRO_XSPEC_VERSION to the "
                "version of XSPEC you are using"
            )

            xspec_version = _default_xspec_version

        else:

            print(f"You have specified the XSPEC version {xspec_version}")

    xspec_version = packaging_version.Version(xspec_version)

    if xspec_version < packaging_version.Version("12.12.0"):
        msg = "WARN: XSPEC version is less than 12.12.0, which is the minimal"
        msg += " supported version for astromodels"
        print(msg)
        return None
    elif xspec_version > packaging_version.Version("12.15.1"):
        msg = "WARN: XSPEC version is greater than 12.15.1, which is the"
        msg += " maximal supported version for astromodels"
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
        (12, 15, 1),
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
    extra_objects = []  # For libraries without unversioned symlinks

    for lib_root in libraries_root:

        this_library, this_library_path, full_lib_path = find_library(
            lib_root, additional_places=[os.path.join(headas_root, "lib")]
        )

        if this_library is None:

            raise IOError(
                "Could not find library %s. Impossible to compile Xspec" % lib_root
            )

        else:

            print("Found library %s in %s" % (this_library, this_library_path))

            if full_lib_path is not None:
                # No unversioned symlink exists, pass the full path directly
                print(
                    "Warning: No unversioned symlink found for %s, "
                    "using full path %s" % (lib_root, full_lib_path)
                )
                extra_objects.append(full_lib_path)
            else:
                # Normal case: unversioned symlink exists, use -l linking
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
    # On macOS, we need to set rpath for libraries in extra_objects
    extra_link_args = []
    if sys.platform.lower().find("darwin") >= 0 and extra_objects:
        # Add rpath for each library directory so macOS can find the libraries
        # at runtime
        for lib_dir in library_dirs:
            extra_link_args.append(f"-Wl,-rpath,{lib_dir}")

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
            extra_link_args=extra_link_args,
            # Add full paths for libraries without unversioned symlinks
            extra_objects=extra_objects,
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
