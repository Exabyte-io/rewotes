from computing.executable import Executable
from computing.b_executable_p2 import BExecutableP2
from qe_error import QEError
import os

class QERun(BExecutableP2):
    """
    Class responsible for actually calling the
    QuantumEspresso executable from within a prepared
    simulation folder.
    
    Various "events" can happen during the execution
    represented by QEError.
    
    Override the handleError(...) method for custom
    error handling.
    """

    PROCESSES_PER_NODE = 1

    def __init__(self, dir, spec, input_file="input.txt"):
        self.dir = dir
        self.spec = spec
        self.input_file = input_file
        self.error_state = QEError.NO_ERROR
        self.stop_flag = False
        self.return_flag = False

    def get_resources(self):
        return self.spec
        
    def get_shell_name(self):
        return "bash"
        
    def get_shell_string(self):
        """
        "mpirun -np " + \
        str(self.spec.get_nodes()*QERun.PROCESSES_PER_NODE) + \
        " pw.x -in " + \
        self.input_file
        """
        return "java -jar /Users/david/Desktop/Update.jar ."
                
        
    def get_std_out(self):
        return os.path.join(self.dir, "stdout.txt")
        
    def get_std_err(self):
        return os.path.join(self.dir, "stderr.txt")
        
    def get_working_dir(self):
        return self.dir
        
    def run(self, envr):
        if not os.path.isdir(self.dir):
            self.handle_error(QEError.SIM_DIR_DOES_NOT_EXIST)
        if self.stop_flag:
            return self.return_flag
        if envr.get_nodes() < self.spec.get_nodes() or \
           envr.get_time() != -1 and spec.get_time() != -1 and \
           envr.get_time() < spec.get_time():
            self.handle_error(QEError.INSUFFICIENT_RESOURCES)
        if self.stop_flag:
            return self.return_flag
        result = super(QERun, self).run(envr)
        if not result:
            self.handle_error(QEError.NONZERO_PROCESS_EXIT)
        if self.stop_flag:
            return self.return_flag
        self.error_state = QEError.NO_ERROR
        return True
        
    def handle_error(self, error):
        self.error_state = error
        if error == QEError.SIM_DIR_DOES_NOT_EXIST:
            print("Fatal Error: Sim Directory does not exist.")
            print(self.dir)
            self.stop_flag = True
            self.return_flag = False
            return
        if error == QEError.INSUFFICIENT_RESOURCES:
            print("Fatal Error: Process was not given sufficient resources.")
            print(self.dir)
            self.stop_flag = True
            self.return_flag = False
            return
        if error == QEError.NONZERO_PROCESS_EXIT:
            print("Fatal Error: Process had a non-zero exit value.")
            print(self.dir)
            self.stop_flag = True
            self.return_flag = False
            return
            
        
        
        
    
