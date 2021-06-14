from .blackbody import (Blackbody, ModifiedBlackbody,
                        NonDissipativePhotosphere,
                        NonDissipativePhotosphere_Deep)
from .functions import (DiracDelta, Exponential_cutoff, Log_parabola, Sin,
                        StepFunction, StepFunctionUpper, has_gsl, has_naima)

if has_naima:
    from .functions import Synchrotron

if has_gsl:

    from .functions import Cutoff_powerlaw_flux

from .absorption import PhAbs, TbAbs, WAbs, has_ebltable
from .apec import has_atomdb

if has_ebltable:
    from .absorption import EBLattenuation

from .apec import has_atomdb
from .extinction import Standard_Rv, ZDust
from .polynomials import (Constant, Cubic, Line, Quadratic, Quartic,
                          get_polynomial)
from .powerlaws import (Band, Band_Calderone, Band_grbm, Broken_powerlaw,
                        Cutoff_powerlaw, Cutoff_powerlaw_Ep,
                        Inverse_cutoff_powerlaw, Powerlaw, Powerlaw_Eflux,
                        Powerlaw_flux, SmoothlyBrokenPowerLaw,
                        Super_cutoff_powerlaw)

if has_atomdb:

    from .apec import APEC, VAPEC



__all__ = ["Band", "Band_Calderone", "Band_grbm", "Broken_powerlaw",
           "Cutoff_powerlaw", "Cutoff_powerlaw_Ep", "Inverse_cutoff_powerlaw", "Powerlaw",
           "Powerlaw_Eflux", "Powerlaw_flux", "SmoothlyBrokenPowerLaw",
           "Super_cutoff_powerlaw",
           "Constant", "Cubic", "DiracDelta", "Exponential_cutoff", "Line",
           "Quadratic", "Sin", "StepFunction", "StepFunctionUpper",
           "PhAbs", "TbAbs", "WAbs",
           "Log_parabola",
           "Blackbody", "ModifiedBlackbody", "NonDissipativePhotosphere", "NonDissipativePhotosphere_Deep",
           "Quartic", "get_polynomial", "ZDust", "Standard_Rv"
           ]

if has_atomdb:
    __all__.extend(["APEC", "VAPEC"])


if has_gsl:

    __all__.extend(["Cutoff_powerlaw_flux"])


if has_naima:

    __all__.extend(["Synchrotron"])


if has_ebltable:

    __all__.extend(["EBLattenuation"])
