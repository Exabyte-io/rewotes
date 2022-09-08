from computing import Executable

class QESetup(Executable):
    """
    Class responsible for preparing a directory for a
    QuantumEspresso simulation, called by QERun.
    
    Various "events" can happen during the execution
    represented by QEError.
    
    Override the handleError(...) method for custom
    error handling.
    """
    
    def get_resources():
        pass
        
    def exec(envr):
        pass
