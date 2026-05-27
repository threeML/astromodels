import math
import os

from astromodels.core.model import Model
from astromodels.core.model_parser import load_model
from astromodels.core.spectral_component import SpectralComponent
from astromodels.core.polarization import (
    LinearPolarization,
    StokesPolarization,
    Unpolarized,
)
from astromodels.functions import Constant, Powerlaw
from astromodels.sources.point_source import PointSource
from pathlib import Path
import numpy as np


def test_linear_polarization_parameters():
    degree = 50.0
    angle = 30.0
    ps = PointSource(
        "PS",
        0,
        0,
        spectral_shape=Powerlaw(),
        polarization=LinearPolarization(degree=degree, angle=angle),
    )
    m1 = Model(ps)
    m1.display()

    m1.save("__test.yml", overwrite=True)

    mp = load_model("__test.yml")
    assert math.isclose(
        mp.sources["PS"].spectrum.main.polarization.degree.value, degree, rel_tol=0.02
    )
    assert math.isclose(
        mp.sources["PS"].spectrum.main.polarization.angle.value, angle, rel_tol=0.02
    )
    mp.display()

    os.remove("__test.yml")


def test_linear_polarization_functions():
    degree = Constant()
    angle = Constant()
    degree.k = 50
    angle.k = 30
    ps = PointSource(
        "PS",
        0,
        0,
        spectral_shape=Powerlaw(),
        polarization=LinearPolarization(degree=degree, angle=angle),
    )
    m1 = Model(ps)
    m1.display()

    m1.save("__test.yml", overwrite=True)

    mp = load_model("__test.yml")
    assert math.isclose(
        mp.sources["PS"].spectrum.main.polarization.degree.Constant.k.value,
        degree.k.value,
        rel_tol=0.02,
    )
    assert math.isclose(
        mp.sources["PS"].spectrum.main.polarization.angle.Constant.k.value,
        angle.k.value,
        rel_tol=0.02,
    )
    mp.display()

    os.remove("__test.yml")


def test_Stokes_polarization_functions():
    u = Constant()
    q = Constant()
    u.k = 0.5
    q.k = 0.5

    ps = PointSource(
        "PS", 0, 0, spectral_shape=Powerlaw(), polarization=StokesPolarization(Q=q, U=u)
    )
    m1 = Model(ps)
    m1.display()

    m1.save("__test.yml", overwrite=True)

    mp = load_model("__test.yml")
    assert math.isclose(
        mp.sources["PS"].spectrum.main.polarization.Q.Constant.k.value,
        q.k.value,
        rel_tol=0.02,
    )
    assert math.isclose(
        mp.sources["PS"].spectrum.main.polarization.Q.Constant.k.value,
        u.k.value,
        rel_tol=0.02,
    )
    mp.display()

    os.remove("__test.yml")

    _ = StokesPolarization(Q=0.5, U=0.5)


def test_unpolarized():
    # should be unpolarized at startupo
    temp_path = Path("__test.yml")
    ps = PointSource("PS", 0, 0, spectral_shape=Powerlaw())

    assert isinstance(
        ps.spectrum.main.polarization, Unpolarized
    ), "Source was not unpolarized after init"
    m1 = Model(ps)
    m1.display()

    m1.save(temp_path, overwrite=True)

    mp = load_model(temp_path)
    assert type(m1.sources["PS"].spectrum.main.polarization) is type(
        mp.sources["PS"].spectrum.main.polarization
    )
    assert isinstance(mp.sources["PS"].spectrum.main.polarization, Unpolarized)

    mp.display()
    temp_path.unlink()

    assert ps.spectrum.main(1, stokes="Q") == ps.spectrum.main(
        1
    ), "Unpolarized changes the value!"


def test_transform():
    u = Constant()
    q = Constant()
    u.k = 0.5
    q.k = 0.5
    polarization = StokesPolarization(Q=q, U=u)
    linear = polarization.to_linear_polarization()
    assert np.isclose(linear(0, stokes="Q"), 0.5)
    assert np.all(
        np.isclose(linear(np.array([0, 1, 2, 3, 4]), stokes="U"), np.ones(5) * 0.5)
    )
    assert np.isclose(linear.degree.value(0), 0.7071067812)


def test_non_callable():
    linear = LinearPolarization(degree=0.2, angle=90)
    spec = SpectralComponent("test", Powerlaw(), polarization=linear)
    spec(0.1, stokes="U")

    assert np.isclose(spec.polarization(0, "U"), 0.0)
    assert np.isclose(spec.polarization(0, "Q"), -0.2)

    stokes = StokesPolarization(Q=0.2, U=0.8)
    spec = SpectralComponent("test", Powerlaw(), polarization=stokes)
    spec(0.1, stokes="U")
    assert spec.polarization(0, "Q") == 0.2
    assert spec.polarization(0, "U") == 0.8
