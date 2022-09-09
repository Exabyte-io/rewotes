from abc import ABC, abstractmethod
import os
import subprocess

class BExecutable(Executable):
    """
    Class automatically handling the details of scheduling a
    SHELL task.
    """
    
    @abstractmethod
    def get_shell_string(self):
        """
        Returns the command that will be interpreted by the SHELL.
        """
        pass
        
    @abstractmethod
    def get_shell_name(self):
        """
        Returns the name of the SHELL that will be invoked to run
        the command.
        """
        pass
        
    @abstractmethod
    def get_std_out(self):
        """
        Returns the absolute canonical filepath for the redirected
        standard out stream. If null, no redirect will occur.
        """
        pass
        
    @abstractmethod
    def get_std_err(self):
        """
        Returns the absolute canonical filepath for the redirected
        standard error stream. If null, no redirect will occur.
        """
        pass
        
    @abstractmethod
    def get_working_dir(self):
        """
        Returns the absolute canonical filepath for the working
        directory in which the command will be run.
        """
        pass
        
    def receive_environment(self, envr):
        """
        Given the current execution envrionment resources,
        sets up the executable to run appropriately.
        """
        pass
        
    def exec(self, envr):
        """
        Runs the specified shell string computational task.
        """
        cmd = self.get_shell_string()
        wd_str = self.get_working_dir()
        if not os.path.isdir(wd_str):
            raise RuntimeException("Unable to locate working directory.")
        std_out = None
        std_err = None
        
        try:
            std_out_call = self.get_std_out()
            if not std_out_call is None:
                std_out = open(std_out_call, "w")
            std_err_call = self.get_std_err()
            if not std_err_call is None:
                std_err = open(std_err_call, "w")
            result = subprocess.run([self.get_shell_name(), "-c", cmd],
                                    stdout=std_out, stderr=std_err,
                                    cwd=wd_str)
            return True
            
        except IOError as e:
            print(e)
            return False
        
        finally:
            if not std_out is None:
                try:
                    std_out.close()
                except IOError as e2:
                    pass
            if not std_err is None:
                try:
                    std_err.close()
                except IOError as e2:
                    pass
        
        
        
    
