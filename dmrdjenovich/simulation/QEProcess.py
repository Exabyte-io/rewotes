from computing import Executable, Resources
import xml.etree.ElementTree as ET
import os

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

    req_resources = Resources(0, 1)

    def __init__(self, dir, spec):
        self.dir = dir
        self.spec = spec
        self.error_state = QEError.NO_ERROR
        self.stop_flag = False
        self.return_flag = False

    def get_dir(self):
        return self.dir
        
    def get_spec(self):
        return self.spec

    def get_resources(self):
        return QEProcess.req_resources
        
    def exec(envr):
        if not os.path.isdir(self.dir):
            self.handle_error(QEError.SIM_DIR_DOES_NOT_EXIST)
        if self.stop_flag:
            return self.return_flag
        
        f_target = os.join(self.dir,
                        self.spec.dict["&CONTROL"].get("prefix", "pwscf") + .xml)
        if not os.path.exists(f_target):
            self.handle_error(QEError.UNABLE_TO_READ_OUTPUT)
        if self.stop_flag:
            return self.return_flag
            
        # Maybe some future code for checking for errors in stdout ...
        
        QEAnalysis analysis = QEAnalysis(f_target)
        get_results(analysis)
        
        self.error_state = QEError.NO_ERROR
        return True
        
    def get_results(self, analysis):
        """
        Method stub meant to be overridden for custom analysis steps.
        """
        return None
        
    def handle_error(self, error):
        self.error_state = error
        if error == QEError.SIM_DIR_DOES_NOT_EXIST:
            print("Fatal Error: Sim Directory not found.")
            print(self.dir)
            self.stop_flag = True
            self.return_flag = False
            return
            
