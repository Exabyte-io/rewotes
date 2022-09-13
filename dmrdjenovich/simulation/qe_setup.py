from computing import Executable, Resources
import os

class QESetup(Executable):
    """
    Class responsible for preparing a directory for a
    QuantumEspresso simulation, called by QERun.
    
    Various "events" can happen during the execution
    represented by QEError.
    
    Override the handleError(...) method for custom
    error handling.
    """
    
    req_resources = Resources(0, -1)
    
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
        return QESetup.req_resources
        
    def run(self, envr):
        if os.path.isdir(self.dir):
            self.handle_error(QEError.DIR_EXISTS)
        if self.stop_flag:
            return self.return_flag
            
        if not os.path.isdir(self.dir):
            try:
                os.path.mkdirs(self.dir)
            except OSError as error:
                self.handle_error(QEError.CANT_CREATE_DIR)
            if self.stop_flag:
                return self.return_flag
        
        spec_path = os.join(self.dir, "input.txt")
        try:
            self.spec.write_to_file(spec_path)
        except IOError:
            self.handle_error(QEError.CANT_WRITE_INPUT_FILES)
        if self.stop_flag:
            return self.return_flag
        self.error_state = QE.NO_ERROR
        return True
        
    def handle_error(self, error):
        self.error_state = error
        if error == QEError.DIR_EXISTS:
            return
        if error == QEError.CANT_CREATE_DIR:
            print("Fatal Error: Can't create the simulation directory.")
            print(self.dir)
            self.stop_flag = True
            self.return_flag = True
            return
        if error = QEError.CANT_WRITE_INPUT_FILES:
            print("Fatal Error: Can't write the input files.")
            print(self.dir)
            self.stop_flag = True
            self.return_flag = True
            return
