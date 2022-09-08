from computing import Executable

class QEProcess(Executable):
    """
    Class that analyzes the outputs of a QuantumEspresso
    simulation.
    
    Pertinent information is extracted from the raw
    outputs and placed in a file called results.txt
    
    Various "events" can happen during the execution
    represented by QEError.
    
    Override the handleError(...) method for custom
    error handling.
    
    Override the getResults(...) method using functions
    from QEAnalysis.py for custom results parsing.
    """

    def get_resources():
        pass
        
    def exec(envr):
        pass
