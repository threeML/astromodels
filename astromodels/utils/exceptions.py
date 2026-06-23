from astromodels.core.my_yaml import my_yaml
from astromodels.utils.disk_usage import disk_usage


class WarningNoTests(ImportWarning):
    pass


class FunctionDefinitionError(Exception):
    pass


class FunctionInstanceError(Exception):
    pass


class DesignViolation(Exception):
    pass


class ModelAssertionViolation(Exception):
    pass


class WrongDimensionality(Exception):
    pass


class TestSpecificationError(Exception):
    pass


class TestFailed(Exception):
    pass


class DocstringIsNotRaw(ValueError):
    pass


class UnknownFunction(ValueError):
    pass


class UnknownParameter(ValueError):
    pass


class EBLTableNotAvailable(ImportWarning):
    pass


class GSLNotAvailable(ImportWarning):
    pass


class NaimaNotAvailable(ImportWarning):
    pass


class InvalidUsageForFunction(Exception):
    pass


class IncompleteGrid(RuntimeError):
    pass


class ValuesNotInGrid(ValueError):
    pass


class MissingDataFile(RuntimeError):
    pass


class InvalidTemplateModelFile(RuntimeError):
    pass


class ModelFileExists(IOError):
    pass


class InvalidInput(ValueError):
    pass


class CannotWriteModel(IOError):
    def __init__(self, directory, message):
        # Add a report on disk usage to the message

        free_space = disk_usage(directory).free

        message += "\nFree space on the file system hosting %s was %.2f Mbytes" % (
            directory,
            free_space / 1024.0 / 1024.0,
        )

        super(CannotWriteModel, self).__init__(message)


class ModelInternalError(ValueError):

    pass


class ModelIOError(IOError):
    pass


class ModelYAMLError(my_yaml.YAMLError):
    pass


class ModelSyntaxError(RuntimeError):
    pass


class WrongCoordinatePair(ValueError):
    pass


class IllegalCoordinateValue(ValueError):
    pass


class WrongCoordinateSystem(ValueError):
    pass


class DuplicatedNode(Exception):
    pass


class ProtectedAttribute(RuntimeError):
    pass


class NonExistingAttribute(RuntimeWarning):
    pass


class UnknownUnit(Exception):
    pass


class UnitMismatch(Exception):
    pass
