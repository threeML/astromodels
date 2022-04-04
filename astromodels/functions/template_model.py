from __future__ import division

import collections
import gc
import os
import re
from builtins import object, range, str
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import astropy.io.fits as fits
import astropy.units as u
import h5py
import numpy as np
import scipy.interpolate
from future.utils import with_metaclass
from interpolation import interp
from interpolation.splines import eval_linear

from astromodels.core.parameter import Parameter
from astromodels.functions.function import Function1D, FunctionMeta
from astromodels.utils import get_user_data_path
from astromodels.utils.logging import setup_logger

log = setup_logger(__name__)

# A very small number which will be substituted to zero during the construction
# of the templates
_TINY_ = 1e-50

__all__ = [
    "IncompleteGrid",
    "ValuesNotInGrid",
    "MissingDataFile",
    "TemplateModelFactory",
    "TemplateModel",
]


class IncompleteGrid(RuntimeError):
    pass


class ValuesNotInGrid(ValueError):
    pass


class MissingDataFile(RuntimeError):
    pass


# This dictionary will keep track of the new classes already created in the current session
_classes_cache = {}


class GridInterpolate(object):
    def __init__(self, grid, values):
        self._grid = grid
        self._values = np.ascontiguousarray(values)

    def __call__(self, v):

        return eval_linear(self._grid, self._values, v)


class UnivariateSpline(object):
    def __init__(self, x, y):

        self._x = x
        self._y = y

    def __call__(self, v):

        return interp(self._x, self._y, v)


class TemplateModelFactory(object):
    def __init__(
        self,
        name: str,
        description: str,
        energies: np.ndarray,
        names_of_parameters: List[str],
        interpolation_degree: int = 1,
        spline_smoothing_factor: int = 0,
    ):

        # Store model name

        # Enforce that it does not contain spaces nor strange characters

        name: str = str(name)

        if re.match("[a-zA-Z_][a-zA-Z0-9_]*", name) is None:
            log.error(
                "The provided name '%s' is not a valid name. You cannot use spaces, "
                "or special characters")
            raise RuntimeError()

        self._name: str = name

        self._description: str = str(description)

        # Store energy grid

        if not isinstance(energies, u.Quantity):

            log.warning(
                "Energy unit is not a Quantity instance, so units has not been provided. Using keV."
            )

            energies: u.Quantity = energies * u.keV

        self._energies: np.ndarray = np.array(energies.to(u.keV).value)

        # Enforce that they are ordered
        self._energies.sort()

        # We create a dictionary which will contain the grid for each parameter

        self._parameters_grids: Dict[
            str, Optional[np.ndarray]] = collections.OrderedDict()

        for parameter_name in names_of_parameters:

            self._parameters_grids[parameter_name] = None

        self._data_frame = None

        self._interpolators = None

        self._interpolation_degree: int = interpolation_degree

        self._spline_smoothing_factor: int = int(spline_smoothing_factor)

    def define_parameter_grid(self, parameter_name: str,
                              grid: np.ndarray) -> None:
        """
        Define the parameter grid for this parameter.
        Pass the name of the parameter and the array of values that it will take in the grid
        """

        if not parameter_name in self._parameters_grids:

            log.error(f"Parameter {parameter_name} is not part of this model")

            raise AssertionError()

        grid_ = np.array(grid)

        if not (grid_.shape[0] > 1):

            log.error(
                "A grid for a parameter must contain at least two elements")

            raise AssertionError()

        # Assert that elements are unique

        if not np.all(np.unique(grid_) == grid_):

            log.error(
                f"Non-unique elements in grid for parameter {parameter_name}")

            raise AssertionError()

        self._parameters_grids[parameter_name] = grid_

    def add_interpolation_data(self, differential_fluxes: np.ndarray,
                               **parameters_values_input: Dict[str, float]):

        # Verify that the grid has been defined for all parameters

        for grid in list(self._parameters_grids.values()):

            if grid is None:

                log.error(
                    "You need to define a grid for all parameters, by using the "
                    "define_parameter_grid method.")

                raise IncompleteGrid()

        # this is the first run and we will create
        # a matrix that stores the data in an
        # n_par,..., n_energies matrix

        if self._data_frame is None:

            shape = []

            for k, v in self._parameters_grids.items():

                shape.append(len(v))

            shape.append(self._energies.shape[0])

            log.debug(f"grid shape: {shape}")

            self._data_frame = np.zeros(tuple(shape))

            log.debug(f"grid shape actual: {self._data_frame.shape}")

            # This is the first data set, create the data frame

            # Create the multi-index

            # self._multi_index = pd.MultiIndex.from_product(
            #     list(self._parameters_grids.values()),
            #     names=list(self._parameters_grids.keys()),
            # )

            # # Pre-fill the data matrix with nans, so we will know if some elements have not been filled

            # self._data_frame = pd.DataFrame(index=self._multi_index,
            #                                 columns=self._energies)

        # Make sure we have all parameters and order the values in the same way as the dictionary

        parameter_idx = []

        for i, (k, v) in enumerate(self._parameters_grids.items()):

            if not k in parameters_values_input:

                log.error(f"Parameter {k} is not in input")

                raise AssertionError()

            parameter_idx.append(
                int(np.where(v == parameters_values_input[k])[0][0]))

        log.debug(f" have index {parameter_idx}")

        # If the user did not specify one of the parameters, then the parameters_values array will contain nan

        if not len(parameter_idx) == len(self._parameters_grids):

            log.error("You didn't specify all parameters' values.")

            raise AssertionError()

        # Make sure we are dealing with pure numpy arrays (list and astropy.Quantity instances will be transformed)
        # First we transform the input into a u.Quantity (if it's not already)

        if not isinstance(differential_fluxes, u.Quantity):

            differential_fluxes = (np.array(differential_fluxes) * 1 /
                                   (u.keV * u.s * u.cm**2))  # type: u.Quantity

        # Then we transform it in the right units and we cast it back to a pure np.array

        differential_fluxes = np.array(
            differential_fluxes.to(1/ (u.keV * u.s * u.cm**2)).value)

        # Now let's check for valid inputs

        if not self._energies.shape[0] == differential_fluxes.shape[0]:

            log.error("Differential fluxes and energies must have "
                      "the same number of elements")

            raise AssertionError()

        # Check that the provided value does not contains nan, inf nor zero (as the interpolation happens in the
        # log space)
        if not np.all(np.isfinite(differential_fluxes)):

            log.error(
                "You have invalid values in the differential flux (nan or inf)"
            )

            raise AssertionError()

        if not np.all(differential_fluxes >= 0):

            log.error(
                "You have negative values in the differential flux (which is of "
                "course impossible)")

            raise AssertionError()

        if not np.all(differential_fluxes > 0):

            log.warning(
                "You have zeros in the differential flux. Since the interpolation happens in the log space, "
                "this cannot be accepted. We will substitute zeros with %g" %
                _TINY_)

            idx = differential_fluxes == 0  # type: np.ndarray
            differential_fluxes[idx] = _TINY_

        # Now set the corresponding values in the data frame

        # Now set the values in the data frame

        self._data_frame[tuple(parameter_idx)] = np.atleast_2d(
            differential_fluxes)

    def save_data(self, overwrite: bool = False):

        # First make sure that the whole data matrix has been filled

        if np.any(np.isnan(self._data_frame)):

            log.error("You have NaNs in the data matrix. Usually this means "
                      "that you didn't fill it up completely, or that some of "
                      "your data contains nans. Cannot save the file.")

            raise AssertionError()

            raise AssertionError()
            
        # Get the data directory

        data_dir_path: Path = get_user_data_path()

        # Sanitize the data file

        filename_sanitized: Path = data_dir_path.absolute() /  f"{self._name}.h5"
        

        # Check that it does not exists
        if filename_sanitized.exists():

            if overwrite:

                try:

                    os.remove(filename_sanitized)

                except:

                    log.error(
                        "The file %s already exists and cannot be removed (maybe you do not have "
                        "permissions to do so?). " % filename_sanitized)

                    raise IOError()

            else:

                log.error(
                    "The file %s already exists! You cannot call two different "
                    "template models with the same name" % filename_sanitized)

                raise IOError()

        # Open the HDF5 file and write objects

        template_file: TemplateFile = TemplateFile(
            name=self._name,
            description=self._description,
            spline_smoothing_factor=self._spline_smoothing_factor,
            interpolation_degree=self._interpolation_degree,
            grid=self._data_frame,
            energies=self._energies,
            parameters=self._parameters_grids,
            parameter_order=list(self._parameters_grids.keys()))

        template_file.save(filename_sanitized)

   

# This adds a method to a class at runtime


def add_method(self, method, name=None):

    if name is None:

        name = method.__name__

    setattr(self.__class__, name, method)


class RectBivariateSplineWrapper(object):
    """
    Wrapper around RectBivariateSpline, which supplies a __call__ method which accept the same
    syntax as the other interpolation methods

    """
    def __init__(self, *args, **kwargs):

        # We can use interp2, which features spline interpolation instead of linear interpolation

        self._interpolator = scipy.interpolate.RectBivariateSpline(
            *args, **kwargs)

    def __call__(self, x):

        res = self._interpolator(*x)

        return res[0][0]


@dataclass
class TemplateFile:
    """
    simple container to read and write 
    the data to an hdf5 file

    """
    name: str
    description: str
    grid: np.ndarray
    parameters: Dict[str, np.ndarray]
    parameter_order: List[str]
    energies: np.ndarray
    interpolation_degree: int
    spline_smoothing_factor: float


    def save(self, file_name: str):

        """
        serialize the contents to a file

        :param file_name: 
        :type file_name: str
        :returns: 

        """
        with h5py.File(file_name, "w") as f:

            f.attrs["name"] = self.name
            f.attrs["description"] = self.description
            f.attrs["interpolation_degree"] = self.interpolation_degree
            f.attrs["spline_smoothing_factor"] = self.spline_smoothing_factor

            f.create_dataset("energies", data=self.energies, compression="gzip")

            f.create_dataset("grid", data=self.grid, compression="gzip")

            # store the parameter order
            dt = h5py.special_dtype(vlen=str)
            po = np.array(self.parameter_order, dtype=dt) 
            f.create_dataset('parameter_order', data=po)
            par_group = f.create_group("parameters")
            for k in self.parameter_order:

                par_group.create_dataset(k, data=self.parameters[k], compression="gzip")

    @classmethod
    def from_file(cls, file_name: str):

        """
        read contents from a file

        :param cls: 
        :type cls: 
        :param file_name: 
        :type file_name: str
        :returns: 

        """
        
        with h5py.File(file_name, "r") as f:

            name = f.attrs["name"]
            description = f.attrs["description"] 
            interpolation_degree = f.attrs["interpolation_degree"] 
            spline_smoothing_factor = f.attrs["spline_smoothing_factor"] 

            energies = f["energies"][()]
            parameter_order = f["parameter_order"][()]

            grid = f["grid"][()]

            parameters = collections.OrderedDict()

            for k in parameter_order:

                parameters[k] = f["parameters"][k][()]
            


        return cls(name=name,
                   description=description,
                   interpolation_degree=interpolation_degree,
                   spline_smoothing_factor=spline_smoothing_factor,
                   energies=energies,
                   parameter_order=parameter_order,
                   parameters = parameters,
                   grid=grid
                   )

                
        


    
class TemplateModel(with_metaclass(FunctionMeta, Function1D)):

    r"""
    description :
        A template model
    latex : $n.a.$
    parameters :
        K :
            desc : Normalization (freeze this to 1 if the template provides the normalization by itself)
            initial value : 1.0
        scale :
            desc : Scale for the independent variable. The templates are handled as if they contains the fluxes
                   at E = scale * x.This is useful for example when the template describe energies in the rest
                   frame, at which point the scale describe the transformation between rest frame energy and
                   observer frame energy. Fix this to 1 to neutralize its effect.
            initial value : 1.0
            min : 1e-5
    """

    def _custom_init_(self, model_name: str, other_name:Optional[str]=None, log_interp: bool=True):
        """
        Custom initialization for this model

        :param model_name: the name of the model, corresponding to the root of the .h5 file in the data directory
        :param other_name: (optional) the name to be used as name of the model when used in astromodels. If None
        (default), use the same name as model_name
        :return: none
        """

        # Get the data directory

        data_dir_path: Path = get_user_data_path()

        # Sanitize the data file

        filename_sanitized = data_dir_path.absolute() /  f"{model_name}.h5"
        

        if not filename_sanitized.exists():

            log.error(f"The data file {filename_sanitized} does not exists. Did you use the "
                "TemplateFactory?")
            
            raise MissingDataFile()

        # Open the template definition and read from it

        self._data_file: Path = filename_sanitized

        # use the file shadow to read
        
        template_file: TemplateFile = TemplateFile.from_file(filename_sanitized)

        self._parameters_grids = collections.OrderedDict()

        

        for key in template_file.parameter_order:

            try:

                # sometimes this is
                # stored binary
                
                k = key.decode()

            except(AttributeError):

                # if not, then we
                # load as a normal str
                
                k = key
                                
            log.debug(f"reading parameter {str(k)}")

  
            self._parameters_grids[str(k)] = template_file.parameters[key]

            

        self._energies = template_file.energies

        # Now get the metadata

        description = template_file.description
        name = template_file.name

        self._interpolation_degree = template_file.interpolation_degree

        self._spline_smoothing_factor = template_file.spline_smoothing_factor

        # Make the dictionary of parameters

        function_definition = collections.OrderedDict()

        function_definition["description"] = description

        function_definition["latex"] = "n.a."

        
        # Now build the parameters according to the content of the parameter grid

        parameters = collections.OrderedDict()

        parameters["K"] = Parameter("K", 1.0)
        parameters["scale"] = Parameter("scale", 1.0)

        for parameter_name in list(self._parameters_grids.keys()):

            grid = self._parameters_grids[parameter_name]

            parameters[parameter_name] = Parameter(
                parameter_name,
                np.median(grid),
                min_value=grid.min(),
                max_value=grid.max(),
            )

        if other_name is None:

            super(TemplateModel, self).__init__(name, function_definition, parameters)

        else:

            super(TemplateModel, self).__init__(
                other_name, function_definition, parameters
            )

        # Finally prepare the interpolators

        self._prepare_interpolators(log_interp, template_file.grid)

        # try can clean up the file
        
        del template_file

        gc.collect()
                
    def _prepare_interpolators(self, log_interp, data_frame):

        # Figure out the shape of the data matrices
        data_shape = [x.shape[0] for x in list(self._parameters_grids.values())]

        self._interpolators = []

        for i, energy in enumerate(self._energies):

            # Make interpolator for this energy
            # NOTE: we interpolate on the logarithm
            # unless specified

            if log_interp:

                this_data = np.array(
                    np.log10(data_frame[...,i]).reshape(*data_shape),
                    dtype=float,
                )

                self._is_log10 = True

            else:

                # work in linear space
                this_data = np.array(
                    data_frame[..., i].reshape(*data_shape), dtype=float
                )

                self._is_log10 = False

                
            if len(list(self._parameters_grids.values())) == 2:

                x, y = list(self._parameters_grids.values())

                # Make sure that the requested polynomial degree is less than the number of data sets in
                # both directions

                msg = (
                    "You cannot use an interpolation degree of %s if you don't provide at least %s points "
                    "in the %s direction. Increase the number of templates or decrease the interpolation "
                    "degree."
                )

                if len(x) <= self._interpolation_degree:

                    log.error(
                        msg
                        % (
                            self._interpolation_degree,
                            self._interpolation_degree + 1,
                            "x",
                        )
                    )

                    raise RuntimeError()

                if len(y) <= self._interpolation_degree:

                    log.error(
                        msg
                        % (
                            self._interpolation_degree,
                            self._interpolation_degree + 1,
                            "y",
                        )
                    )

                    raise RuntimeError()

                this_interpolator = RectBivariateSplineWrapper(
                    x,
                    y,
                    this_data,
                    kx=self._interpolation_degree,
                    ky=self._interpolation_degree,
                    s=self._spline_smoothing_factor,
                )

            else:

                # In more than 2d we can only use linear interpolation

                this_interpolator = GridInterpolate(
                    tuple([np.array(x,dtype='<f8') for x in list(self._parameters_grids.values())]),
                    this_data,
                )

            self._interpolators.append(this_interpolator)

        # clear the data
        self._data_frame = None
            
    def _set_units(self, x_unit, y_unit):

        self.K.unit = y_unit

        self.scale.unit = 1 / x_unit

    # This function will be substituted during construction by another version with
    # all the parameters of this template

    def evaluate(self, x, K, scale, *args):

        return K * self._interpolate(x, scale, args)

    def _interpolate(self, energies, scale, parameters_values):

        if isinstance(energies, u.Quantity):

            # Templates are always saved with energy in keV. We need to transform it to
            # a dimensionless quantity (actually we take the .value property) because otherwise
            # the logarithm below will fail.

            energies = np.array(
                energies.to("keV").value, ndmin=1, copy=False, dtype=float
            )

            # Same for the scale

            scale = scale.to(1/ u.keV).value

        if self._is_log10:

            log_energies = np.log10(energies)

        else:

            log_energies = energies

        e_tilde = self._energies * scale

        # Gather all interpolations for these parameters' values at all defined energies
        # (these are the logarithm of the values)
        # note that if these are not logged, then the name is superflous

        log_interpolations = np.array(
            [
                self._interpolators[i](np.atleast_1d(parameters_values))
                for i in range(self._energies.shape[0])
            ]
        )

        # Now interpolate the interpolations to get the flux at the requested energies

        # NOTE: the variable "interpolations" contains already the log10 of the values,

        if self._is_log10:

            interpolator = UnivariateSpline(np.log10(e_tilde), log_interpolations)

            values = np.power(10., interpolator(log_energies))

        else:

            interpolator = UnivariateSpline(e_tilde, log_interpolations)

            values = interpolator(log_energies)

        # The division by scale results from the differential:
        # E = e * scale
        # de = dE / scale
        # dN / dE = dN / de * de / dE = dN / de * (1 / scale)

        # NOTE: the units are added back through the multiplication by K in the evaluate method

        return values / scale

    def clean(self):
        """
        Table models can consume a lot of memory. If are creating lots of 
        table models in memory for simulations, you may want to call
        clean on the model try and remove some of the memory consumed by the models

        :returns: 

        """
        
        self._interpolators = None
        del self._interpolators
        gc.collect()

        log.info("You have 'cleaned' the table model and it will no longer be useable")

    # def __del__(self):

    #     self.clean()

        
    @property
    def data_file(self):

        return self._data_file

    def to_dict(self, minimal=False):

        data = super(Function1D, self).to_dict(minimal)

        # if not minimal:
        #
        #     data['extra_setup'] = {'data_file': self._data_file}

        return data


class XSPECTableModel(object):
    def __init__(
        self,
        xspec_table_model_file,
        interpolation_degree=1,
        spline_smoothing_factor=0,
        log_centers=True,
    ):
        """
        Convert an XSPEC table model to an astromodels TemplateModel.
        usage:

        xtm = XSPECTableModel("ST95.fits")
        xtm.to_table_model('test', 'test')

        reloaded_table_model = TemplateModel('test', log_interp=False)

        Note: if the reloaded model is returning NaNs, adjust the interpolation
        scheme

        :param xspec_table_model_file:
        :param interpolation_degree:
        :param spline_smoothing_factor: spline smoothing
        :param log_centers: treat energies with log centers
        :returns:
        :rtype:

        """

        self._interpolation_degree = interpolation_degree
        self._spline_smoothing_factor = spline_smoothing_factor
        self._log_centers = log_centers

        self._xspec_file_name = xspec_table_model_file
        self._extract_model()

    def _extract_model(self):

        with fits.open(self._xspec_file_name) as f:

            # get the energies

            energies = f["ENERGIES"]
            ene_lo = energies.data["ENERG_LO"]
            ene_hi = energies.data["ENERG_HI"]
            spectra = f["SPECTRA"]
            
            if f[0].header["MODLUNIT"] == 'photons/cm^2/s' :
               log.info("This table model has flux units of photons/cm^2/s. It will be converted to differential Flux by dividing by the energy bin widths. Logarithmic centers will be turned off. ")
               delta_ene = ene_hi-ene_lo
               self._spectrum = spectra.data["INTPSPEC"]/delta_ene
               self._log_centers = False
            else:
               self._spectrum = spectra.data["INTPSPEC"]
            
            # log centers

            if self._log_centers:
                self._energy = np.sqrt(ene_lo * ene_hi)

            else:

                self._energy = (ene_hi + ene_lo) / 2.0

            params = f["PARAMETERS"]

            self._names = params.data["NAME"]

            self._n_params = len(self._names)

            self._params_dict = {}

            for i, name in enumerate(self._names):

                this_dict = {}

                this_dict["pmin"] = params.data["MINIMUM"][i]
                this_dict["pmax"] = params.data["MAXIMUM"][i]
                if self._n_params > 1:

                    try:

                        this_dict["values"] = spectra.data["PARAMVAL{%d}" % i]

                    except:
                        this_dict["values"] = spectra.data["PARAMVAL"][:, i]

                else:

                    this_dict["values"] = spectra.data["PARAMVAL"]

                self._params_dict[name] = this_dict

    def to_table_model(self, file_name, model_name, overwrite=False):
        """
        Write the table model to your local astromodels database

        :param file_name: name of file to store
        :param model_name: name of the model
        :param overwrite: overwite the previous model
        :returns:
        :rtype:

        """

        tmf = TemplateModelFactory(
            file_name,
            model_name,
            self._energy,
            self._names,
            self._interpolation_degree,
            self._spline_smoothing_factor,
        )

        for name, param_table in self._params_dict.items():

            tmf.define_parameter_grid(name, np.unique(param_table["values"]))

        for i in range(self._spectrum.shape[0]):

            input_dict = {}
            for k, v in self._params_dict.items():
                input_dict[k] = v["values"][i]

            tmf.add_interpolation_data(self._spectrum[i, :], **input_dict)

        tmf.save_data(overwrite=overwrite)



def convert_old_table_model(model_name: str):
    from pandas import HDFStore

    # Get the data directory

    data_dir_path: Path = get_user_data_path()

    # Sanitize the data file
    
    filename_sanitized = data_dir_path / f"{model_name}.h5"
    
    
    if not filename_sanitized.exists():

        log.error("The data file %s does not exists. Did you use the "
             "TemplateFactory?" % (filename_sanitized))
         
        raise MissingDataFile( )
    
    with HDFStore(filename_sanitized) as store:
    
        data_frame = store["data_frame"]

        parameters_grids = collections.OrderedDict()

        processed_parameters = 0

        for key in list(store.keys()):

            match = re.search("p_([0-9]+)_(.+)", key)
            
            if match is None:
                
                continue
            
            else:
                
                tokens = match.groups()
                
                this_parameter_number = int(tokens[0])
                this_parameter_name = str(tokens[1])
                
                
                if not this_parameter_number == processed_parameters:

                    log.error("Parameters out of order!")

                    raise AssertionError()
                
                parameters_grids[this_parameter_name] = store[key]
                
                processed_parameters += 1
                
        energies = np.array(store["energies"])

        shape = [len(x) for _, x in parameters_grids.items()]

        shape.append(len(energies))

        grid = np.zeros(shape)

        for i, e in enumerate(energies):

            data = data_frame[e].values.reshape(*shape[:-1])

            grid[...,i] = data
        
                     
        
        # Now get the metadata
            
        metadata = store.get_storer("data_frame").attrs.metadata
        
        description = metadata["description"]
        
        name = metadata["name"]
        
        interpolation_degree = metadata["interpolation_degree"]
        
        spline_smoothing_factor = metadata["spline_smoothing_factor"]
        
        # Make the dictionary of parameters

        function_definition = collections.OrderedDict()

        function_definition["description"] = description



    template_file: TemplateFile = TemplateFile(
        name=name,
        description=description,
        spline_smoothing_factor=spline_smoothing_factor,
        interpolation_degree=interpolation_degree,
        grid=grid,
        energies=energies,
        parameters=parameters_grids,
        parameter_order=list(parameters_grids.keys()))

    template_file.save(filename_sanitized)


        
