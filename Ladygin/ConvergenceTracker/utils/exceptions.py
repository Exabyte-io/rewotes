class OutdirInconsistencyWarning(Warning):
    """Warning for output dir is not working dir"""
    pass

class NofileWarning(Warning):
    """Warning for absence of input setting file"""
    pass

class NotConvergedWarning(Warning):
    """Warning for kpoint search not being converged in the range of available kpoints"""
    pass

class MissingPseudoError(Exception):
    """Missing pseudo error"""
    pass