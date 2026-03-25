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
