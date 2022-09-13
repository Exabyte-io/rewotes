from enum import Enum

class QEError(Enum):
    """
    Enum determining the error-state of an actual QuantumEspresso
    simulation.
    """
    
    NO_ERROR = 1
    
    DIR_EXISTS = 2
    CANT_CREATE_DIR = 3
    CANT_WRITE_INPUT_FILES = 4
    
    SIM_DIR_DOES_NOT_EXIST = 5
    INSUFFICIENT_RESOURCES = 6
    NONZERO_PROCESS_EXIT = 7
    
    UNABLE_TO_READ_OUTPUT = 8
    
    
    
