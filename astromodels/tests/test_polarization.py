import math
import os

from astromodels.core.model import Model
from astromodels.core.model_parser import load_model
from astromodels.core.polarization import (
    LinearPolarization,
    StokesPolarization,
    Unpolarized,
)
from astromodels.functions import Constant, Powerlaw
from astromodels.sources.point_source import PointSource


def test_linear_polarization_parameters(tmp_path):
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
    temp_file = tmp_path / "__test.yml"
    m1.save(temp_file, overwrite=True)

    mp = load_model(temp_file)
    assert math.isclose(
        mp.sources["PS"].spectrum.main.polarization.degree.value, degree, rel_tol=0.02
    )
    assert math.isclose(
        mp.sources["PS"].spectrum.main.polarization.angle.value, angle, rel_tol=0.02
    )
    mp.display()

    os.remove(temp_file)


def test_linear_polarization_functions(tmp_path):
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
    temp_file = tmp_path / "__test.yml"
    m1.save(temp_file, overwrite=True)

    mp = load_model(temp_file)
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

    os.remove(temp_file)


def test_Stokes_polarization_functions(tmp_path):
    u = Constant()
    q = Constant()
    u.k = 0.5
    q.k = 0.5

    ps = PointSource(
        "PS", 0, 0, spectral_shape=Powerlaw(), polarization=StokesPolarization(Q=q, U=u)
    )
    m1 = Model(ps)
    m1.display()
    temp_file = tmp_path / "__test.yml"
    m1.save(temp_file, overwrite=True)

    mp = load_model(temp_file)
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

    os.remove(temp_file)

    _ = StokesPolarization(Q=0.5, U=0.5)


def test_unpolarized(tmp_path):
    # should be unpolarized at startupo
    temp_path = tmp_path / "__test.yml"
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
