from __future__ import division

import os

import copy

import numpy as np
import math

from astromodels import u, update_logging_level
from astromodels.core.polarization import *
from astromodels.core.model import Model
from astromodels.core.model_parser import *
from astromodels.functions import Constant,Powerlaw
from astromodels.sources.point_source import PointSource

update_logging_level("DEBUG")

def test_linear_polarization_parameters():
    degree = 50.0
    angle = 30.0
    ps = PointSource('PS',0,0,spectral_shape=Powerlaw(),polarization=LinearPolarization(degree=degree, angle=angle))
    m1 = Model(ps)
    m1.display()

    m1.save("__test.yml",overwrite=True)

    mp = load_model("__test.yml")
    assert math.isclose(mp.sources["PS"].spectrum.main.polarization.degree.value, degree, rel_tol=0.02)
    assert math.isclose(mp.sources["PS"].spectrum.main.polarization.angle.value, angle, rel_tol=0.02)
    mp.display()

    os.remove("__test.yml")

def test_linear_polarization_functions():
    degree = Constant()
    angle  = Constant()
    degree.k = 50
    angle.k  = 30
    ps = PointSource('PS',0,0,spectral_shape=Powerlaw(),polarization=LinearPolarization(degree=degree, angle=angle))
    m1 = Model(ps)
    m1.display()

    m1.save("__test.yml",overwrite=True)

    mp = load_model("__test.yml")
    assert math.isclose(mp.sources["PS"].spectrum.main.polarization.degree.Constant.k.value, degree.k.value, rel_tol=0.02)
    assert math.isclose(mp.sources["PS"].spectrum.main.polarization.angle.Constant.k.value, angle.k.value, rel_tol=0.02)
    mp.display()

    os.remove("__test.yml")

def test_Stokes_polarization_functions():
    u = Constant()
    q = Constant()
    u.k = 0.5
    q.k = 0.5

    ps = PointSource('PS',0,0,spectral_shape=Powerlaw(),polarization=StokesPolarization(Q=q, U=u))
    m1 = Model(ps)
    m1.display()

    m1.save("__test.yml",overwrite=True)

    mp = load_model("__test.yml")
    assert math.isclose(mp.sources["PS"].spectrum.main.polarization.Q.Constant.k.value, q.k.value, rel_tol=0.02)
    assert math.isclose(mp.sources["PS"].spectrum.main.polarization.Q.Constant.k.value, u.k.value, rel_tol=0.02)
    mp.display()

    os.remove("__test.yml")