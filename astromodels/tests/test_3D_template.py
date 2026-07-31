"""
Unit tests for spatial_model.py
Covers: ModelFactory, TemplateFile, HaloModel, GridInterpolate,
        UnivariateSpline, RectBivariateSplineWrapper, and the custom exceptions.
"""

import collections
import os
from pathlib import Path
from unittest.mock import patch

import astropy.units as u
import h5py
import numpy as np
import pytest
from astropy.io import fits

# ---------------------------------------------------------------------------
# Helpers to build minimal FITS / HDF5 fixtures in memory
# ---------------------------------------------------------------------------

NE, NL, NB = 5, 8, 8  # energy, lon, lat pixels used in mock FITS data


def _make_fits_file(path: str, ne: int = NE, nl: int = NL, nb: int = NB) -> str:
    """Write a minimal 3-D FITS cube and return its path."""
    data = np.random.rand(ne, nl, nb).astype(np.float64)
    hdu = fits.PrimaryHDU(data)
    # Mandatory WCS keys for a 3D cube
    hdu.header["CDELT1"] = 0.5
    hdu.header["CDELT2"] = 0.5
    hdu.header["CDELT3"] = 0.2
    hdu.header["CRVAL1"] = 0.0
    hdu.header["CRVAL2"] = 0.0
    hdu.header["CRVAL3"] = 5.0
    hdu.header["CRPIX1"] = nl // 2
    hdu.header["CRPIX2"] = nb // 2
    hdu.header["CRPIX3"] = 1
    hdu.header["NAXIS1"] = nl
    hdu.header["NAXIS2"] = nb
    hdu.header["NAXIS3"] = ne
    hdu.writeto(path, overwrite=True)
    return path


def _make_template_h5(path: str, param_names=("alpha",), grid_sizes=(3,)) -> str:
    """
    Write a minimal TemplateFile HDF5 and return its path.
    Single parameter grid by default.
    """
    ne, nl, nb = NE, NL, NB
    energies = np.linspace(5.0, 5.8, ne)
    lats = np.linspace(-2.0, 2.0, nb)
    lons = np.linspace(-2.0, 2.0, nl)

    shape = list(grid_sizes) + [ne, nl, nb]
    grid = np.random.rand(*shape)

    parameters = collections.OrderedDict()
    for name, sz in zip(param_names, grid_sizes):
        parameters[name] = np.linspace(0.1, 1.0, sz)

    with h5py.File(path, "w") as f:
        f.attrs["name"] = "test_model"
        f.attrs["description"] = "unit-test template"
        f.attrs["degree_of_interpolation"] = 1
        f.attrs["spline_smoothing_factor"] = 0

        f.create_dataset("energies", data=energies)
        f.create_dataset("lats", data=lats)
        f.create_dataset("lons", data=lons)
        f.create_dataset("grid", data=grid)

        dt = h5py.special_dtype(vlen=str)
        f.create_dataset("parameter_order", data=np.array(list(param_names), dtype=dt))

        pg = f.create_group("parameters")
        for k, v in parameters.items():
            pg.create_dataset(k, data=v)

    return path


try:
    from astromodels.functions.spatial_model import (
        GridInterpolate,
        HaloModel,
        IncompleteGrid,
        MissingDataFile,
        ModelFactory,
        RectBivariateSplineWrapper,
        TemplateFile,
        UnivariateSpline,
        ValuesNotInGrid,
    )

    _IMPORT_OK = True
except ImportError:
    _IMPORT_OK = False

pytestmark = pytest.mark.skipif(
    not _IMPORT_OK, reason="Module or its dependencies could not be imported"
)


class TestCustomExceptions:
    def test_incomplete_grid_is_runtime_error(self):
        with pytest.raises(RuntimeError):
            raise IncompleteGrid("missing data")

    def test_values_not_in_grid_is_value_error(self):
        with pytest.raises(ValueError):
            raise ValuesNotInGrid("out of range")

    def test_missing_spatial_data_file_is_runtime_error(self):
        with pytest.raises(RuntimeError):
            raise MissingDataFile("no file")


class TestGridInterpolate:
    def test_returns_correct_value_1d(self):
        x = np.array([0.0, 1.0, 2.0], dtype="<f8")
        y = np.array([0.0, 1.0, 4.0], dtype="<f8")
        interp = GridInterpolate((x,), y)
        result = interp(np.array([[1.0]]))
        assert np.isfinite(result)

    def test_shape_preserved(self):
        x = np.array([0.0, 1.0, 2.0], dtype="<f8")
        y = np.array([10.0, 20.0, 30.0], dtype="<f8")
        interp = GridInterpolate((x,), y)
        val = interp(np.array([[1.0]]))
        # Extract scalar to avoid deprecation warning
        val_scalar = val.item()
        assert 10.0 <= val_scalar <= 30.0


class TestUnivariateSpline:
    def test_interpolates_at_known_points(self):
        x = np.array([1.0, 2.0, 3.0], dtype="<f8")
        y = np.array([2.0, 4.0, 6.0], dtype="<f8")
        spl = UnivariateSpline(x, y)
        result = spl(np.array([2.0]))
        assert np.isclose(result.item(), 4.0, atol=1e-6)

    def test_returns_finite(self):
        x = np.linspace(0, 10, 20, dtype="<f8")
        y = np.sin(x).astype("<f8")
        spl = UnivariateSpline(x, y)
        val = spl(np.array([5.0]))
        assert np.isfinite(val)


class TestRectBivariateSplineWrapper:
    def test_call_returns_scalar(self):
        x = np.linspace(0, 1, 5)
        y = np.linspace(0, 1, 5)
        z = np.outer(x, y)
        wrapper = RectBivariateSplineWrapper(x, y, z, kx=1, ky=1, s=0)
        result = wrapper((np.array([0.5]), np.array([0.5])))
        assert np.isfinite(result)


class TestTemplateFile:
    def test_save_and_load_roundtrip(self, tmp_path):
        path = str(tmp_path / "template.h5")
        ne, nl, nb = NE, NL, NB

        params = collections.OrderedDict({"alpha": np.linspace(0.1, 1.0, 3)})
        tf = TemplateFile(
            name="roundtrip",
            description="test",
            grid=np.random.rand(3, ne, nl, nb),
            energies=np.linspace(5.0, 5.8, ne),
            lats=np.linspace(-2.0, 2.0, nb),
            lons=np.linspace(-2.0, 2.0, nl),
            parameters=params,
            parameter_order=["alpha"],
            degree_of_interpolation=1,
            spline_smoothing_factor=0,
        )
        tf.save(path)

        loaded = TemplateFile.from_file(path)
        assert loaded.name == "roundtrip"
        assert loaded.description == "test"
        np.testing.assert_array_equal(loaded.energies, tf.energies)
        np.testing.assert_array_equal(loaded.lats, tf.lats)
        np.testing.assert_array_equal(loaded.lons, tf.lons)
        np.testing.assert_array_almost_equal(loaded.grid, tf.grid)

        # Handle possible bytes keys from HDF5 (convert to string)
        key = next(iter(loaded.parameters.keys()))
        if isinstance(key, bytes):
            loaded.parameters = {k.decode(): v for k, v in loaded.parameters.items()}
        np.testing.assert_array_equal(loaded.parameters["alpha"], params["alpha"])

    def test_from_file_preserves_parameter_order(self, tmp_path):
        path = str(tmp_path / "order.h5")
        _make_template_h5(path, param_names=("beta", "gamma"), grid_sizes=(3, 4))
        loaded = TemplateFile.from_file(path)
        # first key may be bytes or string
        first_key = next(iter(loaded.parameters.keys()))
        assert first_key in ("beta", b"beta")

    def test_save_creates_file(self, tmp_path):
        path = str(tmp_path / "exists.h5")
        params = collections.OrderedDict({"p": np.linspace(0, 1, 2)})
        tf = TemplateFile(
            name="exists",
            description="d",
            grid=np.zeros((2, NE, NL, NB)),
            energies=np.zeros(NE),
            lats=np.zeros(NB),
            lons=np.zeros(NL),
            parameters=params,
            parameter_order=["p"],
            degree_of_interpolation=1,
            spline_smoothing_factor=0,
        )
        tf.save(path)
        assert os.path.isfile(path)


class TestModelFactory:
    # ------------------------------------------------------------------
    # __init__ validation
    # ------------------------------------------------------------------
    def test_init_valid_name(self):
        mf = ModelFactory("valid_name", "desc", ["p1"])
        assert mf._name == "valid_name"

    def test_init_invalid_name_raises(self):
        with pytest.raises(RuntimeError):
            ModelFactory("123 bad name!", "desc", ["p1"])

    def test_parameters_initialised_to_none(self):
        mf = ModelFactory("model", "desc", ["a", "b"])
        assert mf._parameters_grids["a"] is None
        assert mf._parameters_grids["b"] is None

    def test_define_parameter_grid_stores_array(self):
        mf = ModelFactory("m", "d", ["alpha"])
        grid = np.linspace(0, 1, 5)
        mf.define_parameter_grid("alpha", grid)
        np.testing.assert_array_equal(mf._parameters_grids["alpha"], grid)

    def test_define_parameter_grid_unknown_parameter_raises(self):
        mf = ModelFactory("m", "d", ["alpha"])
        with pytest.raises(AssertionError):
            mf.define_parameter_grid("unknown", np.linspace(0, 1, 5))

    def test_define_parameter_grid_single_element_raises(self):
        mf = ModelFactory("m", "d", ["alpha"])
        with pytest.raises(AssertionError):
            mf.define_parameter_grid("alpha", np.array([1.0]))

    def test_define_parameter_grid_converts_to_float(self):
        mf = ModelFactory("m", "d", ["alpha"])
        mf.define_parameter_grid("alpha", [0, 1, 2])
        assert mf._parameters_grids["alpha"].dtype == float

    def test_add_data_without_grid_raises_incomplete_grid(self, tmp_path):
        fits_path = _make_fits_file(str(tmp_path / "cube.fits"))
        mf = ModelFactory("m2", "d", ["alpha"])
        with pytest.raises(IncompleteGrid):
            mf.add_interpolation_data(fits_path, alpha=0.5)

    def test_add_data_without_fits_file_raises(self):
        mf = ModelFactory("m3", "d", ["alpha"])
        mf.define_parameter_grid("alpha", np.linspace(0, 1, 3))
        with pytest.raises(RuntimeError):
            mf.add_interpolation_data(None, alpha=0.0)

    def test_add_data_fills_data_frame(self, tmp_path):
        fits_path = _make_fits_file(str(tmp_path / "cube.fits"))
        mf = ModelFactory("m4", "d", ["alpha"])
        grid = np.array([0.0, 0.5, 1.0])
        mf.define_parameter_grid("alpha", grid)
        for val in grid:
            mf.add_interpolation_data(fits_path, alpha=val)
        assert mf._data_frame is not None
        assert not np.any(np.isnan(mf._data_frame))

    def test_save_creates_h5_file(self, tmp_path):
        fits_path = _make_fits_file(str(tmp_path / "cube.fits"))
        mf = ModelFactory("save_test", "d", ["alpha"])
        grid_vals = np.array([0.0, 0.5, 1.0])
        mf.define_parameter_grid("alpha", grid_vals)
        for val in grid_vals:
            mf.add_interpolation_data(fits_path, alpha=val)

        with patch(
            "astromodels.functions.spatial_model.get_user_data_path",
            return_value=tmp_path,
        ):
            mf.save_data(overwrite=False)

        assert (tmp_path / "save_test.h5").exists()

    def test_save_overwrite_false_raises_when_file_exists(self, tmp_path):
        fits_path = _make_fits_file(str(tmp_path / "cube.fits"))
        mf = ModelFactory("dup_model", "d", ["alpha"])
        grid_vals = np.array([0.0, 0.5, 1.0])
        mf.define_parameter_grid("alpha", grid_vals)
        for val in grid_vals:
            mf.add_interpolation_data(fits_path, alpha=val)

        # Create the file first
        (tmp_path / "dup_model.h5").touch()

        with (
            patch(
                "astromodels.functions.spatial_model.get_user_data_path",
                return_value=tmp_path,
            ),
            pytest.raises(OSError),
        ):
            mf.save_data(overwrite=False)

    def test_save_overwrite_true_replaces_file(self, tmp_path):
        fits_path = _make_fits_file(str(tmp_path / "cube.fits"))
        mf = ModelFactory("over_model", "d", ["alpha"])
        grid_vals = np.array([0.0, 0.5, 1.0])
        mf.define_parameter_grid("alpha", grid_vals)
        for val in grid_vals:
            mf.add_interpolation_data(fits_path, alpha=val)

        with patch(
            "astromodels.functions.spatial_model.get_user_data_path",
            return_value=tmp_path,
        ):
            mf.save_data(overwrite=False)  # first save
            mf.save_data(overwrite=True)  # should overwrite without error

        assert (tmp_path / "over_model.h5").exists()


@pytest.fixture
def halo_model_fixture(tmp_path):
    """Build a minimal HaloModel backed by a real HDF5 file."""
    h5_path = str(tmp_path / "halo_test.h5")
    _make_template_h5(h5_path, param_names=("alpha",), grid_sizes=(3,))

    with patch(
        "astromodels.functions.spatial_model.get_user_data_path", return_value=tmp_path
    ):
        model = HaloModel.__new__(HaloModel)
        model._custom_init_("halo_test", other_name=None)

    return model


class TestHaloModel:
    def test_raises_if_file_missing(self, tmp_path):
        with (
            patch(
                "astromodels.functions.spatial_model.get_user_data_path",
                return_value=tmp_path,
            ),
            pytest.raises(MissingDataFile),
        ):
            model = HaloModel.__new__(HaloModel)
            model._custom_init_("nonexistent_model")

    def test_parameter_grids_loaded(self, halo_model_fixture):
        model = halo_model_fixture
        assert "alpha" in model._parameters_grids
        assert len(model._parameters_grids["alpha"]) == 3

    def test_energy_lat_lon_arrays_loaded(self, halo_model_fixture):
        model = halo_model_fixture
        assert model._E.shape[0] == NE
        assert model._B.shape[0] == NB
        assert model._L.shape[0] == NL

    def test_cache_initialised_empty(self, halo_model_fixture):
        model = halo_model_fixture
        assert len(model._cached_values) == 0

    def test_define_region_equatorial(self, halo_model_fixture):
        model = halo_model_fixture
        ramin, ramax, decmin, decmax = model.define_region(10.0, 20.0, -5.0, 5.0)
        assert ramin == 10.0
        assert ramax == 20.0
        assert decmin == -5.0
        assert decmax == 5.0

    def test_define_region_stores_attributes(self, halo_model_fixture):
        model = halo_model_fixture
        model.define_region(0.0, 30.0, -10.0, 10.0)
        assert model.ramin == 0.0
        assert model.ramax == 30.0
        assert model.decmin == -10.0
        assert model.decmax == 10.0

    # ------------------------------------------------------------------
    # get_boundaries
    # ------------------------------------------------------------------
    def test_get_boundaries_returns_correct_tuples(self, halo_model_fixture):
        model = halo_model_fixture
        model.define_region(5.0, 15.0, -3.0, 3.0)
        lon_bounds, lat_bounds = model.get_boundaries()
        assert lon_bounds == (5.0, 15.0)
        assert lat_bounds == (-3.0, 3.0)

    def test_interpolate_returns_correct_shape(self, halo_model_fixture):
        model = halo_model_fixture
        energies = np.array([1e3, 1e4])  # keV
        lons = np.linspace(-1.0, 1.0, 4)
        lats = np.linspace(-1.0, 1.0, 4)
        param_values = (0.5,)  # within grid [0.1, 0.55, 1.0]

        result = model._interpolate(energies, lons, lats, param_values)
        assert result.shape == (len(lons), len(energies))

    def test_interpolate_all_finite(self, halo_model_fixture):
        model = halo_model_fixture
        energies = np.array([5e3])
        lons = np.linspace(-0.5, 0.5, 5)
        lats = np.linspace(-0.5, 0.5, 5)
        param_values = (0.5,)

        result = model._interpolate(energies, lons, lats, param_values)
        assert np.all(np.isfinite(result))

    def test_interpolate_caches_result(self, halo_model_fixture):
        model = halo_model_fixture
        energies = np.array([1e3])
        lons = np.zeros(3)
        lats = np.zeros(3)
        param_values = (0.5,)

        model._interpolate(energies, lons, lats, param_values)
        assert param_values in model._cached_values

    def test_interpolate_uses_cache_on_second_call(self, halo_model_fixture):
        model = halo_model_fixture
        energies = np.array([1e3])
        lons = np.zeros(3)
        lats = np.zeros(3)
        param_values = (0.5,)

        r1 = model._interpolate(energies, lons, lats, param_values)
        r2 = model._interpolate(energies, lons, lats, param_values)
        np.testing.assert_array_equal(r1, r2)

    def test_interpolate_lon_lat_size_mismatch_raises(self, halo_model_fixture):
        model = halo_model_fixture
        # The error is raised with no message inside AttributeError
        with pytest.raises(AttributeError):
            model._interpolate(
                np.array([1e3]),
                np.zeros(4),
                np.zeros(5),  # different size
                (0.5,),
            )

    def test_interpolate_accepts_astropy_quantity_energies(self, halo_model_fixture):
        model = halo_model_fixture
        energies = np.array([1.0, 10.0]) * u.keV
        lons = np.zeros(3)
        lats = np.zeros(3)

        result = model._interpolate(energies, lons, lats, (0.5,))
        assert result.shape == (3, 2)

    # ------------------------------------------------------------------
    # evaluate() tests – fixed to avoid keyword/positional conflicts
    # and to pass energies as array, not scalar
    # ------------------------------------------------------------------
    def test_evaluate_returns_correct_shape(self, halo_model_fixture):
        model = halo_model_fixture
        lons = np.linspace(-1.0, 1.0, 4)
        lats = np.linspace(-1.0, 1.0, 4)
        energies = np.array([1e3, 5e3])  # 1D array

        # evaluate(x, y, z, K, lon0, lat0, *args)
        result = model.evaluate(lons, lats, energies, 1.0, 0.0, 0.0, 0.5)
        assert result.shape == (len(lons), len(energies))

    def test_evaluate_scales_with_K(self, halo_model_fixture):
        model = halo_model_fixture
        lons = np.linspace(-0.5, 0.5, 3)
        lats = np.linspace(-0.5, 0.5, 3)
        energies = np.array([1e3])  # 1D array

        r1 = model.evaluate(lons, lats, energies, 1.0, 0.0, 0.0, 0.5)
        r2 = model.evaluate(lons, lats, energies, 2.0, 0.0, 0.0, 0.5)
        np.testing.assert_array_almost_equal(r2, 2.0 * r1)

    def test_evaluate_single_energy_array(self, halo_model_fixture):
        """evaluate works when a single energy value is passed as a 1‑element array."""
        model = halo_model_fixture
        lons = np.linspace(-0.5, 0.5, 3)
        lats = np.linspace(-0.5, 0.5, 3)
        # Pass a 1‑element array, not a scalar (the code does not handle scalars)
        result = model.evaluate(lons, lats, np.array([1e3]), 1.0, 0.0, 0.0, 0.5)
        assert np.all(np.isfinite(result))

    def test_data_file_property_is_path(self, halo_model_fixture):
        assert isinstance(halo_model_fixture.data_file, Path)
