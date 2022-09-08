from computing import Executable

class QERun(Executable):
    """
    Class responsible for actually calling the
    QuantumEspresso executable from within a prepared
    simulation folder.
    
    Various "events" can happen during the execution
    represented by QEError.
    
    Override the handleError(...) method for custom
    error handling.
    """

    def get_resources():
        pass
        
    def exec(envr):
        pass
