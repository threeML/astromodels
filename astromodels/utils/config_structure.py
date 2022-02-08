import logging
from dataclasses import dataclass
from enum import IntEnum, Enum


# from .catalog_structure import Catalogs, PublicDataServer
# from .fitting_structure import BayesianDefault, MLEDefault
# from .plotting_structure import GenericPlotting, ModelPlotting
# from .plugin_structure import Plugins, TimeSeries
# from .point_source_structure import PointSourceDefaults

# logging
class LoggingLevel(IntEnum):
    DEBUG = logging.DEBUG
    INFO = logging.INFO
    WARNING = logging.WARNING
    ERROR = logging.ERROR
    CRITICAL = logging.CRITICAL


@dataclass
class Logging:

    path: str = "~/.astromodels/log"
    developer: bool = 'off'
    usr: bool = 'on'
    console: bool = 'on'
    level: LoggingLevel = LoggingLevel.INFO
    startup_warnings: bool = 'on'


class AbsTables(Enum):
    WILM = "WILM"
    ASPL = "ASPL"
    AG89 = "AG89" 


class EBLTable(Enum):
    franceschini = "franceschini"
    kneiske = "kneiske"
    dominguez = "dominguez"
    inuoe = "inuoe"
    gilmore = "gilmore"

    
    
@dataclass
class AbsorptionModels:
    tbabs_table: AbsTables = AbsTables.WILM
    phabs_table: AbsTables = AbsTables.AG89
    ebl_table: EBLTable = EBLTable.dominguez

@dataclass
class Modeling:
    use_memoization: bool = True
    use_parameter_transforms: bool = True
    ignore_parameter_bounds: bool = False




@dataclass
class Config:
    logging: Logging = Logging()
    absorption_models: AbsorptionModels = AbsorptionModels()
    modeling: Modeling = Modeling()
