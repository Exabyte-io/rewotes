"""Custom exceptions."""


class InputFilesNotFoundError(Exception):
    """Error when the provided input path string is not a valid Path."""


class TargetOutputNotFoundError(Exception):
    """Error when requested target variable was not found in solver results."""


class SolverSubprocessFailedError(Exception):
    """Return processed subprocess error."""


class SolverNotImplementedError(Exception):
    """Error when the requested solver is not implemented."""


class SolverNotInstalledError(Exception):
    """Error when the requested solver is not installed."""


class UnsupportedParameterError(Exception):
    """Error when provided parameters are unsupported."""
