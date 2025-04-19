class OutdirInconsistencyWarning(Warning):
    """Warning for output dir is not working dir"""
    pass

class NofileWarning(Warning):
    """Warning for absence of input setting file"""
    pass

class MissingPseudoError(Exception):
    """Missing pseudo error"""
    pass