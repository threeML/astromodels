"""Template fitting of 3D energy dependent model"""

import collections
import gc
import hashlib
import os
import re
from builtins import range, str
from dataclasses import dataclass
from itertools import product
from pathlib import Path

import astropy.units as u
import h5py
import numpy as np
import scipy.interpolate
from astropy.coordinates import SkyCoord
from astropy.io import fits
from numpy.typing import NDArray
from scipy.interpolate import RectBivariateSpline, RegularGridInterpolator
from typing_extensions import List, Optional, OrderedDict, Sequence, Tuple, TypeAlias

from astromodels.core.parameter import Parameter
from astromodels.functions.function import Function3D, FunctionMeta
from astromodels.utils import get_user_data_path
from astromodels.utils.logging import setup_logger

ndarray: TypeAlias = NDArray[np.float64]

log = setup_logger(__name__)

__author__ = "Ramiro"
__version__ = "0.1"
__comment__ = "Aims to be part of the astromodels package"


# NOTE: Script adapted GalProp and TemplateModelFactory in Astromodels.

__all__ = [
    "IncompleteGrid",
    "ValuesNotInGrid",
    "MissingSpatialDataFile",
    "ModelFactory",
    "HaloModel",
]


class IncompleteGrid(RuntimeError):
    """Check that the grid has been correctly filled in the add_data() method

    :raises if: Raised if the grid contains any Nans or None values
    :type RuntimeError: Exception
    """


class ValuesNotInGrid(ValueError):
    """Check whether there are values not contained within the user defined grid

    :raises if: There are values not contained within the grid
    :type ValueError: Exception
    """


class MissingSpatialDataFile(RuntimeError):
    """Checks if the file exists

    :raises if: File is not present
    :type RuntimeError: Exception
    """


# This dictionary will keep track of the new classes already
# created in the current session
_classes_cache = {}


class GridInterpolate:
    """Interpolation over a regular grid of n dimension (limited to linear
    interpolation)"""

    def __init__(self, grid: ndarray, values: ndarray) -> None:
        self._grid = grid
        self._values = np.ascontiguousarray(values)

    def __call__(self, v) -> None:
        return scipy.interpolate.interpn(self._grid, self._values, v)


class UnivariateSpline:
    """Simple one dimensional spline interpolation"""

    def __init__(self, x, y) -> None:
        self._x = x
        self._y = y

    def __call__(self, v):
        if isinstance(self._x, np.ndarray) and isinstance(self._y, np.ndarray):
            return np.interp(v, self._x, self._y.reshape(self._x.shape))
        else:
            np.interp(v, self._x, self._y)


# This class builds a dataframe from morphology parameters for a source
# The information from the source comes from FITS files with
# 3D template model maps
class ModelFactory:
    def __init__(
        self,
        name: str,
        description: str,
        names_of_parameters: List[str],
        degree_of_interpolation: int = 1,
        spline_smoothing_factor: int = 0,
    ) -> None:
        """Class builds a data table from 3d mapcube templates with information
        on morphology and spectral parameters for different energy bins. These
        parameters are used for interpolation in the class HaloModel

        :param name: Name of output HDF5 file
        :param description: A brief summary of the template model for reference
        :param names_of_parameters: unique name for parameters to interpolate over
        :param degree_of_interpolation: Polynomial degree for interpolating,
        defaults to 1
        :param spline_smoothing_factor: Smoothing factor used during spline
        interpolation, defaults to 0
        :raises RuntimeError: Raised if name of the file cannot contain spaces or
        special characters.
        Example:
        >>> model_factory = ModelFactory(
        ...     "my_model",
        ...     "A brief summary of the model",
        ...     ["param1", "param2"],
        ...     degree_of_interpolation=1,
        ...     spline_smoothing_factor=0,
        ... )
        >>> model_factory.define_parameter_grid(
        ...     "param1", np.linspace(0, 10, 100)
        ... )
        >>> model_factory.define_parameter_grid(
        ...     "param2", np.linspace(0, 10, 100)
        ... )
        >>> for i, val1 in enumerate(model_factory._parameters_grids["param1"]):
        >>>     for j, val2 in enumerate(model_factory._parameters_grids["param2"]):
        >>>         model_factory.add_interpolation_data(fits_file="template.fits",
        ...                                             param1=val1,
        ...                                             param2=val2)
        >>> # Save the data to a file.
        ... # If the data file already exists, it won't overwrite
        ... # the current file until overwrite is set to overwrite = True
        >>> model_factory.save_data()
        """

        self._data_frame: Optional[ndarray] = None
        self._delLon: float
        self._delLat: float
        self._delEn: float
        self._refLon: float
        self._refLat: float
        self._refEn: float
        self._refLonPix: float
        self._refLatPix: float
        self._refEnPix: float
        self._nl: int
        self._nb: int
        self._ne: int
        self._map: ndarray
        self._L: ndarray
        self._B: ndarray
        self._E: ndarray

        # Store the model name
        # Ensure that it contains no spaces nor special characters
        name = str(name)

        if re.match("[a-zA-Z_][a-zA-Z0-9_]*", name) is None:
            log.error(
                f"The provided name {name} is not a valid name. "
                "You cannot use spaces, nor special characters."
            )

            raise RuntimeError()

        self._name = name
        self._description = str(description)
        self._degree_of_interpolation = int(degree_of_interpolation)
        self._spline_smoothing_factor = int(spline_smoothing_factor)

        # Create a dictionary which will contain the grid for each parameter
        self._parameters_grids: OrderedDict[str, ndarray] = collections.OrderedDict()

        for parameter_name in names_of_parameters:
            self._parameters_grids[parameter_name] = None  # type: ignore

    def define_parameter_grid(self, parameter_name: str, grid: ndarray) -> None:
        """Defines the user provider parameter grid for a given parameter with its
        associated name

        :param parameter_name: Name of parameter for later access when using model
        :param grid: Array of allowed values for the paremter
        :raises AssertionError: Check that the parameter is part of the model
        :raises AssertionError: Ensure all parameter values in the grid are
        non-repeating
        """

        if parameter_name not in self._parameters_grids:
            log.error(f"Parameter {parameter_name} is not part of this model.")

            raise AssertionError()

        # if the grid is not numpy array, conver it to one
        grid_ = np.array(grid, dtype=float)

        if grid_.shape[0] <= 1:
            log.error(
                "A grid for a parameter must containt at least two "
                "elements for interpolation."
            )

            raise AssertionError()

        # Assert that all elements are unique
        if not np.all(np.unique(grid_) == grid_):
            log.error(
                f"Non-unique elements found in grid of parameter {parameter_name}."
            )

        self._parameters_grids[parameter_name] = grid_

    def add_interpolation_data(
        self, fits_file: str, ihdu: int = 0, **parameters_values_input: float
    ) -> None:
        """Fill data table with information from 3D mapcube templates.

        :param fits_file: FITS file with 3D mapcube
        :param ihdu: Path to primary HDU, defaults to 0
        :raises IncompleteGrid: Check if the grid is not filled with strange values
        :raises RuntimeError: Checks that a FITS file was provided
        :raises AssertionError: Ensure we have the right number of parameters as
        provided in the define_parameter_grid method
        :return: none
        """

        # Verify that a grid has been defined for all parameters
        for grid in list(self._parameters_grids.values()):
            if grid is None:
                log.error(
                    "You need to define a grid for all parameters, "
                    "by using the define_parameter_grid method."
                )

                raise IncompleteGrid()

        # Now load the information from the FITS 3D template model
        # template models contain flux values in terms of energy, RA, and dec

        if fits_file is None:
            log.error("You need to specify a FITS file with a template map.")

            raise RuntimeError()

        self._fits_file: Path = Path(fits_file)

        with fits.open(self._fits_file) as f:
            self._delLon = f[ihdu].header["CDELT1"]  # type: ignore
            self._delLat = f[ihdu].header["CDELT2"]  # type: ignore
            self._delEn = 0.2  # f[ihdu].header["CDELT3"]
            self._refLon = f[ihdu].header["CRVAL1"]  # type: ignore
            self._refLat = f[ihdu].header["CRVAL2"]  # type: ignore
            self._refEn = 5  # f[ihdu].header["CRVAL3"]  # Log(E/MeV) -> GeV to MeV
            self._refLonPix = f[ihdu].header["CRPIX1"]  # type: ignore
            self._refLatPix = f[ihdu].header["CRPIX2"]  # type: ignore
            self._refEnPix = f[ihdu].header["CRPIX3"]  # type: ignore

            # 3D array containing the flux values
            self._map: ndarray = np.array(f[ihdu].data, dtype=float)  # type: ignore

            self._nl = f[ihdu].header["NAXIS1"]  # Longitude # type: ignore
            self._nb = f[ihdu].header["NAXIS2"]  # Latitude # type: ignore
            self._ne = f[ihdu].header["NAXIS3"]  # Energy # type: ignore

            self._L = np.linspace(
                self._refLon - self._refLonPix * self._delLon,
                self._refLon + (self._nl - self._refLonPix - 1) * self._delLon,
                self._nl,
                dtype=float,
            )

            self._B = np.linspace(
                self._refLat - self._refLatPix * self._delLat,
                self._refLat + (self._nb - self._refLatPix - 1) * self._delLat,
                self._nb,
                dtype=float,
            )

            self._E = np.linspace(
                self._refEn,
                self._refEn + (self._ne - 1) * self._delEn,
                self._ne,
                dtype=float,
            )

            # for i, e in enumerate(self._E):
            #
            #    self._map[i] = np.fliplr(self._map[i])

            h = hashlib.sha224()
            h.update(self._map)
            self.hash = int(h.hexdigest(), 16)

        if self._data_frame is None:
            shape = [len(v) for k, v in self._parameters_grids.items()]

            shape.extend([self._E.shape[0], self._L.shape[0], self._B.shape[0]])

            log.debug(f"grid shape: {shape}")

            self._data_frame = np.zeros(tuple(shape))

            log.debug(f"grid shape actual: {self._data_frame.shape}")

            # This is the first data set, create the data frame
            # The dataframe is indexed by a parameter multi-index for rows and
            # with the 3D template multi-index information as the columns

            # log.info("Creating the multi-indices....")

        # Make sure we have all parameters and order the values in the same way
        # as the dictionary
        parameter_idx = []
        for l, (key, val) in enumerate(self._parameters_grids.items()):
            if key not in parameters_values_input:
                log.error(f"Parameter {key} is not in input")

            parameter_idx.append(
                int(np.where(val == parameters_values_input[key])[0][0])
            )

        log.debug(f"have index {parameter_idx}")

        if len(parameter_idx) != len(self._parameters_grids):
            log.error("You didn't specify all parameters' values")

            raise AssertionError()

        # filling the data map
        self._data_frame[tuple(parameter_idx)] = self._map

    def save_data(self, overwrite: bool = False) -> None:
        """Save the table into a file for later usage with SpatialModel.

        :param overwrite: Allows the overwriting of an already existing file,
        defaults to False
        :raises AssertionError: Raised if there are any strange values within the table
        :raises IOError: Raised if the filed exists and deleting it not permitted
        :raises IOError: Raised if the filed exists and overriting it is not enabled
        """

        # First make sure that the whole data matrix has been filled
        if np.any(np.isnan(self._data_frame)):
            log.error(
                "You have NaNs in the data matrix. Usually this means "
                "that you didn't fill it up completely, or that some of "
                "your data contains nans. Cannot save the file."
            )

            raise ValueError()

        # Get the data directory
        data_dir_path: Path = get_user_data_path()

        # Sanitize the data file
        filename_sanitized: Path = data_dir_path / f"{self._name}.h5"

        # Check if the file already exists
        # if os.path.exists(filename_sanitized):
        if filename_sanitized.exists():
            if overwrite:
                try:
                    os.remove(filename_sanitized)

                except Exception as exc:
                    log.error(
                        f"The file {filename_sanitized} already exists. "
                        "and cannot be removed (maybe you do not have "
                        "enough permissions to do so?)."
                    )

                    raise IOError() from exc

            else:
                log.error(
                    f"The file {filename_sanitized} already exists! "
                    "You cannot call two different spatial models with the "
                    "same name."
                )

                raise IOError()

        if self._data_frame is not None:
            # Open the HDF5 and write objects
            template_file: TemplateFile = TemplateFile(
                name=self._name,
                description=self._description,
                spline_smoothing_factor=self._spline_smoothing_factor,
                degree_of_interpolation=self._degree_of_interpolation,
                grid=self._data_frame,
                energies=self._E,
                lats=self._B,
                lons=self._L,
                parameters=self._parameters_grids,
                parameter_order=list(self._parameters_grids.keys()),
            )

            template_file.save(filename_sanitized.as_posix())
        else:
            raise RuntimeError("No data frame to save")


# This adds a method to a class at run time
def add_method(self, method, name=None) -> None:
    """Add a method to a class at run time

    :param method: function to add to the class
    :param name: method name, defaults to None
    """
    if name is None:
        name = method.func_name

    setattr(self.__class__, name, method)


class RectBivariateSplineWrapper:
    """Wrapper class around RectBivariateSpline which supplies a __call__ method
    which accepts the same syntax as other interpolation methods
    """

    def __init__(self, *args, **kwargs) -> None:
        # We can use interp2, which features spline interpolation instead of
        # linear interpolation
        self._interpolator = scipy.interpolate.RectBivariateSpline(*args, **kwargs)

    def __call__(self, x):
        res = self._interpolator(*x)

        return res[0][0]


@dataclass
class TemplateFile:
    """
    simple container to read and write the
    data to an hdf5 file

    """

    name: str
    description: str
    grid: ndarray
    energies: ndarray
    lats: ndarray
    lons: ndarray
    parameters: OrderedDict[str, ndarray]
    parameter_order: List[str]
    degree_of_interpolation: int
    spline_smoothing_factor: int

    def save(self, file_name: str) -> None:
        """Serialize the contents to a file and save it

        :param file_name: Path of output file
        :return: none
        """

        with h5py.File(file_name, "w") as f:
            f.attrs["name"] = self.name
            f.attrs["description"] = self.description
            f.attrs["degree_of_interpolation"] = self.degree_of_interpolation
            f.attrs["spline_smoothing_factor"] = self.spline_smoothing_factor

            f.create_dataset("energies", data=self.energies, compression="gzip")
            f.create_dataset("lats", data=self.lats, compression="gzip")
            f.create_dataset("lons", data=self.lons, compression="gzip")
            f.create_dataset("grid", data=self.grid, compression="gzip")

            # store the parameter order
            dt = h5py.special_dtype(vlen=str)
            po = np.array(self.parameter_order, dtype=dt)
            f.create_dataset("parameter_order", data=po)
            par_group = f.create_group("parameters")
            for k in self.parameter_order:
                par_group.create_dataset(k, data=self.parameters[k], compression="gzip")

    @classmethod
    def from_file(cls, file_name: str):
        """Read the contents from a template file

        :param file_name: Path to the file
        :return: TemplateFile holding the data from the file
        """
        with h5py.File(file_name, "r") as f:
            name: str = f.attrs["name"]
            description: str = f.attrs["description"]
            degree_of_interpolation: int = f.attrs["degree_of_interpolation"]
            spline_smoothing_factor: int = f.attrs["spline_smoothing_factor"]

            parameter_order: List[str] = f["parameter_order"][()]  # type: ignore
            energies: ndarray = f["energies"][()]  # type: ignore
            lats: ndarray = f["lats"][()]  # type: ignore
            lons: ndarray = f["lons"][()]  # type: ignore

            grid: ndarray = f["grid"][()]  # type: ignore

            parameters: OrderedDict[str, ndarray] = collections.OrderedDict()

            for k in parameter_order:
                parameters[k] = f["parameters"][k][()]  # type: ignore

        return cls(
            name=name,
            description=description,
            degree_of_interpolation=degree_of_interpolation,
            spline_smoothing_factor=spline_smoothing_factor,
            energies=energies,
            lats=lats,
            lons=lons,
            parameter_order=parameter_order,
            parameters=parameters,
            grid=grid,
        )


# class SpatialModel(with_metaclass(FunctionMeta, Function3D)):
class HaloModel(Function3D, metaclass=FunctionMeta):
    r"""
    description: 3D interpolation over morphology of a source using FITS templates
                with spectral and spatial information.

    latex: $n.a.$


    parameters:

        K:

            desc: Normalization (freeze this to 1 if the template provides the
                    normalization by itself)
            initial value: 1
            fix: yes

        lon0:

            desc: Longitude of the center of source
            initial value: 0.0
            min: 0.0
            max: 360.0

        lat0:

            desc: Latitude of the center of source
            initial value: 0.0
            min: -90.0
            max: 90.0
    """

    def _custom_init_(self, model_name: str, other_name: Optional[str] = None) -> None:
        """
        Custom initialization for this model
        :param model_name: the name of the model, corresponding to the
        root of the .h5 file in the data directory
        :param other_name: (optional) the name to be used as name of the model
        when used in astromodels. If None
        (default), use the same as model_name.
        :return: none
        """

        # Get the data directory
        data_dir_path: Path = get_user_data_path()

        # Sanitize the file
        filename_sanitized = data_dir_path.absolute() / f"{model_name}.h5"

        if not filename_sanitized.exists():
            log.error(
                f"The data file {filename_sanitized} does not exist."
                " Did you use the ModelFactory?"
            )

            raise MissingSpatialDataFile

        # Open the template definition and read from it
        self._data_file = filename_sanitized

        # use the file shadow to read
        template_file: TemplateFile = TemplateFile.from_file(
            filename_sanitized.as_posix()
        )

        self._parameters_grids = collections.OrderedDict()

        for key in template_file.parameter_order:
            try:
                # sometimes this is stored binary
                k = key.decode()

            except AttributeError:
                # if not, load as a normal str
                k = key

            log.debug(f"reading parameter {str(k)}")

            self._parameters_grids[str(k)] = template_file.parameters[key]

        # get the template parameters
        self._E: ndarray = template_file.energies
        self._B: ndarray = template_file.lats
        self._L: ndarray = template_file.lons

        # get the dataframe
        # self.grid = template_file.grid

        # Now get the metadata
        description = template_file.description
        name = template_file.name

        self._degree_of_interpolation = template_file.degree_of_interpolation
        self._spline_smoothing_factor = template_file.spline_smoothing_factor

        # Make the dictionary of parameters for the model
        function_definition: collections.OrderedDict[str, str] = (
            collections.OrderedDict()
        )
        function_definition["description"] = description
        function_definition["latex"] = "n.a."

        # Now build the parameters according to the content of the parameter grid
        parameters: collections.OrderedDict[str, Parameter] = collections.OrderedDict()

        parameters["K"] = Parameter("K", 1.0)
        parameters["lon0"] = Parameter("lon0", 0.0, min_value=0.0, max_value=360.0)
        parameters["lat0"] = Parameter("lat0", 0.0, min_value=-90.0, max_value=90.0)

        for parameter_name in list(self._parameters_grids.keys()):
            grid = self._parameters_grids[parameter_name]

            parameters[parameter_name] = Parameter(
                parameter_name,
                np.median(grid),
                min_value=grid.min(),
                max_value=grid.max(),
            )

        if other_name is None:
            super().__init__(name, function_definition, parameters)

        else:
            super().__init__(other_name, function_definition, parameters)

        self._setup()

        self._data_frame = template_file.grid

        # self._prepare_interpolators(log_interp=True, data_frame=template_file.grid)
        # clean things up a bit

        # Setup cache to avoid unnecessary computations
        self._cached_values: OrderedDict[Tuple[float, ...], ndarray] = (
            collections.OrderedDict()
        )

        del template_file

        gc.collect()

    @staticmethod
    def _msg(interp_degree: int, parameter_name: str) -> str:
        """Retrun a message for the user

        :param interp_degree: degree of interpolation
        :param parameter_name: name of parameter over which to interpolate
        :return: message to user
        """
        return (
            f"You cannot use an interpolation degree of {interp_degree} if "
            f"you don't provide at least {interp_degree} points "
            f"in the {parameter_name} direction. Increase the number of "
            "templates or decrease interpolation degree."
        )

    def _prepare_interpolators(self, log_interp: bool, data_frame: ndarray) -> None:
        """Reads data frame of normalized flux values and performs interpolation
        over morphology parameters declared in ModelFactory's parameter grid

        :param log_interp: interpolation carried on log scale
        :param data_frame: Data table with information from template grid
        :raises RuntimeError: Raised if degree of interpolation for x is not at
        least > 1 than number of parameters in x direction
        least > 1 than number of parameters in y direction
        :return: none
        """

        log.info("Preparing the interpolators...")

        # Figure out the shape of the data matrices
        para_shape = np.array(
            [x.shape[0] for x in list(self._parameters_grids.values())]
        )
        parameter_grid_values: List[float] = list(self._parameters_grids.values())
        parameter_values_len: int = len(parameter_grid_values)

        # interpolate over the parameters
        self._interpolators: List[Optional[RectBivariateSpline, GridInterpolate]] = []

        this_interpolator: Optional[
            UnivariateSpline, RectBivariateSpline, GridInterpolate
        ] = None

        indices = product(
            range(self._E.shape[0]), range(self._L.shape[0]), range(self._B.shape[0])
        )

        for i, j, k in indices:
            reshaped_data = np.array(
                data_frame[..., i, j, k].reshape(*para_shape), dtype=float
            )
            if log_interp:
                reshaped_data = np.log10(reshaped_data)
                self._is_log10 = True
            else:
                self._is_log10 = False

            this_data = reshaped_data

            if parameter_values_len == 1:
                xpoints = np.array(
                    [parameter_grid_values[x] for x in range(parameter_values_len)][0]
                )
                ypoints = np.array([this_data[x] for x in range(this_data.shape[0])])

                this_interpolator = UnivariateSpline(xpoints, ypoints)

            elif parameter_values_len == 2:
                x, y = list(self._parameters_grids.values())

                # Make sure that the requested polynomial degree is
                # less thant the number of data sets in
                # both directions

                if len(x) <= self._degree_of_interpolation:
                    log.error(self._msg(self._degree_of_interpolation, "x"))

                    raise RuntimeError()

                if len(y) <= self._degree_of_interpolation:
                    log.error(self._msg(self._degree_of_interpolation, "y"))

                    raise RuntimeError()

                this_interpolator = RectBivariateSplineWrapper(  # type: ignore
                    x,
                    y,
                    this_data,
                    kx=self._degree_of_interpolation,
                    ky=self._degree_of_interpolation,
                    s=self._spline_smoothing_factor,
                )

            else:
                # In more than 2d, we can only interpolate linearly
                this_interpolator = GridInterpolate(
                    tuple(
                        [np.array(x, dtype="<f8") for x in self._parameter_grid_values]
                    ),
                    this_data,
                )

            self._interpolators.append(this_interpolator)

        del data_frame
        gc.collect()

    def _interpolate(
        self,
        energies: ndarray,
        lons: ndarray,
        lats: ndarray,
        parameter_values: tuple[float, ...],
    ) -> ndarray:
        """Evaluates the morphology parameters and generates interpolated map
        that is then used for the interpolation over energy RA, Dec and energy

        :param energies: Energy values
        :param lons: Longitudes within the extended source boundaries
        :param lats: Latitutdes within the extended source boundaries
        :param parameter_values: User provided morphology parameters defined in
        ModelFactory
        :raises AttributeError: Longitudes and latitudes do not have the same
        dimensions
        :return: Map of interpolated values over energies, longitudes, and latitudes
        """

        if not hasattr(self, "_interpolators"):
            self._prepare_interpolators(log_interp=False, data_frame=self._data_frame)

        if isinstance(energies, u.Quantity):
            energies = np.array(
                energies.to("keV").value, ndmin=1, copy=False, dtype=float
            )

        # transform energy from keV to MeV (galprop likes MeV, 3ML likes keV)
        log_energies = np.log10(energies * (1.0 * u.keV).to(u.MeV).value)

        if self._cached_values.get(parameter_values) is not None:
            # self._log_interpolated_map = self._cached_values[parameter_values]
            return self._cached_values[parameter_values]

        evaluated_map: Sequence = list(
            map(
                lambda interpolator: interpolator(np.atleast_1d(parameter_values)),
                self._interpolators,
            )
        )

        self._log_interpolated_map = np.array(evaluated_map, dtype="<f8")

        # limit the size of the cache to 30 values and pop the oldest entries
        # inserted valued (follows the FIFO approach)
        if len(self._cached_values) > 30:
            while len(self._cached_values) > 10:
                self._cached_values.popitem(last=False)

        map_shape = [x.shape[0] for x in [self._E, self._L, self._B]]

        self._interpolator = RegularGridInterpolator(
            (self._E, self._L, self._B),
            self._log_interpolated_map.reshape(*map_shape),
            bounds_error=False,
        )

        if lons.size != lats.size:
            log.error("Lons and lats should have the same size!")

            raise AttributeError()

        # f_interpolated = np.zeros([energies.size, lats.size])
        f_interpolated = np.zeros([lons.size, energies.size])

        # # evaluate the interpolators over energy, ra, and dec
        for i, e in enumerate(log_energies):
            engs: ndarray = np.repeat(e, lats.shape[0])

            interpolated_slice = self._interpolator((engs, lons, lats))

            # clean up the data
            bad_idx = np.isnan(interpolated_slice)
            interpolated_slice[bad_idx] = 0
            bad_idx = np.isinf(interpolated_slice)
            interpolated_slice[bad_idx] = 0

            if self._is_log10:
                # NOTE: if interpolation is carried using the log10 scale,
                # ensure that values outside range of interpolation remain
                # zero after conversion to linear scale.
                # because if function returns zero, 10**(0) = 1.
                # This affects the fit in 3ML and breaks things.

                log_interpolated_slice = np.array(
                    [0.0 if x == 0.0 else np.power(10.0, x) for x in interpolated_slice]
                )
                f_interpolated[:, i] = log_interpolated_slice

            else:
                f_interpolated[:, i] = interpolated_slice

        assert np.all(np.isfinite(f_interpolated)), "some values are wrong!"

        self._cached_values[parameter_values] = f_interpolated

        return f_interpolated

    def _setup(self):
        self._frame = "icrs"  # ICRS()

    def clean(self) -> None:
        """
        Table models can consume a lof of memory.
        This method calls a clean method and removes some of the memory
        consumed by the models.
        :returns: none
        """
        self._interpolators = None  # type: ignore
        del self._interpolators
        gc.collect()

        log.info("You have cleaned the table model and it will no longer be usable.")

    def __del__(self) -> None:
        # clear the cache
        # self._cached_values.clear()
        self.clean()

    def _set_units(self, x_unit, y_unit, z_unit, w_unit) -> None:
        self.lon0.unit = x_unit
        self.lat0.unit = y_unit

        # self.K.unit = 1 / (u.MeV * u.cm**2 * u.s * u.sr)
        # keep this units to if templates have been normalized
        self.K.unit = 1 / (u.sr)

    def evaluate(self, x, y, z, K, lon0, lat0, *args):
        lons = x
        lats = y
        energies = z

        # angsep = angular_distance_fast(lon0, lat0, lons, lats)

        # if only one energy is passed, make sure we can iterate just once
        if not isinstance(energies, np.ndarray):
            energies = np.array(energies)

        interpolated_image = self._interpolate(energies, lons, lats, args)

        # return np.multiply(K, self._interpolate(log_energies, lons, lats, args))
        # return np.multiply(K, interpolated_image/1000)  # type: ignore
        # to go from MeV to keV
        # return np.multiply(K, interpolated_image / 1000)

        # if templates are normalized no need to convert back
        return np.multiply(K, interpolated_image)  # type: ignore

    # def set_frame(self, new_frame):
    # """
    # Set a new frame for the coordinates (the default is ICRS J2000)
    # :param new_frame: a coordinate frame from astropy
    # :return: (none)
    # """
    #
    # assert isinstance(new_frame, BaseCoordinateFrame)
    #
    # self._frame = new_frame

    @property
    def data_file(self):
        return self._data_file

    def to_dict(self, minimal: bool = False):
        data = super(Function3D, self).to_dict(minimal)

        if not minimal:
            data["extra_setup"] = {
                "_frame": self._frame,
                "ramin": self.ramin,
                "ramax": self.ramax,
                "decmin": self.decmin,
                "decmax": self.decmax,
            }

        return data

    # Define the region within the template ROI
    def define_region(
        self,
        xcoord_start: float,
        xcoord_end: float,
        ycoord_start: float,
        ycoord_end: float,
        galactic: bool = False,
    ) -> tuple[float, float, float, float]:
        """Defined boundaries of template

        :param xcoord_start: Minimum longitude (RA)
        :param xcoord_end: Maximum longitude (RA)
        :param ycoord_start: Minimum latitute (Dec)
        :param ycoord_end: Maximum latitue (Dec)
        :param galactic: Determine whether coordinates are galactic, defaults to False
        :return: coordinates of the boundaries of the template
        """

        if galactic:
            lmin: float = xcoord_start
            lmax: float = xcoord_end
            bmin: float = ycoord_start
            bmax: float = ycoord_end

            _coord = SkyCoord(
                l=[lmin, lmin, lmax, lmax], b=[bmin, bmax, bmax, bmin], frame="galactic"
            )

            self.ramin: float = min(_coord.transform_to(self.frame.value).ra.value)
            self.ramax: float = max(_coord.transform_to(self.frame.value).ra.value)
            self.decmin: float = min(_coord.transform_to(self.frame.value).dec.value)
            self.decmax: float = max(_coord.transform_to(self.frame.value).dec.value)

        else:
            self.ramin = xcoord_start
            self.ramax = xcoord_end
            self.decmin = ycoord_start
            self.decmax = ycoord_end

        return self.ramin, self.ramax, self.decmin, self.decmax

    def get_boundaries(self) -> tuple[tuple[float, float], tuple[float, float]]:
        """Get the boundaries of the template

        :return: Boundaries of template for (min longitude, max longtiude),
        (min latitude, max latitude)
        """
        min_longitude: float = self.ramin
        max_longitude: float = self.ramax
        min_latitude: float = self.decmin
        max_latitude: float = self.decmax

        return (min_longitude, max_longitude), (min_latitude, max_latitude)
