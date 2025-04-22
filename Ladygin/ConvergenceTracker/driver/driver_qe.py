from ase.io import read
from ase.io import espresso
import warnings
import os


from ..utils import NofileWarning, MissingPseudoError, OutdirInconsistencyWarning
from ..utils import pw_template, atoms_template, parse_qe_input, convert_settings
from .calculator import calculator_qe

from .driver_base import Driver


class DriverQE(Driver):
    """Driver for point calculations using quantum espresso"""
    
    def __init__(self, 
                 workdir: str = './',
                 calculator:calculator_qe = calculator_qe(),
                 input_file_name: str = 'pw.in',
                 encut: float = 40,
                 **kwargs) -> None:
        """input_file_name - simulation settings input file name
           calculator - job runner for point calculation
           workdir - working directory
           encut - energy cutoff for point calc
        """
        super().__init__(workdir, calculator, input_file_name, encut, **kwargs)

        self.input_file = os.path.join(workdir, input_file_name)
    
    
    def gen_input(self, kpoint: float) -> None:
        """Generates the input file for the driver based on input parameters
        
            encut - energy cutoff (in Ry)
            kpoint - kpoint dimension    
        """
        if not os.path.exists(self.input_file):
            warnings.warn("The input file is missing I hope you know what you are doing. Calcs would be done silicon reference", NofileWarning)
            atoms = atoms_template
            settings = pw_template
        else:
            atoms = espresso.read_espresso_in(f'{self.workdir}/pw.in')
            settings = parse_qe_input(f'{self.workdir}/pw.in')

        # Convert to ase format
        settings_ase, pseudo = convert_settings(settings)
        
        for key in pseudo:
            if not os.path.exists(os.path.join(self.workdir, "pseudo", pseudo[key])):
                raise MissingPseudoError(f"Missing pseudo {pseudo[key]}, calculations will not run!")


        try:
            assert settings_ase['outdir'] == self.workdir and settings_ase['pseudo_dir'] == os.path.join(self.workdir, "pseudo")
        except AssertionError:
            warnings.warn("Output dir or pseudo_dir is inconsistant with working dir. Replacing outdir with working dir", OutdirInconsistencyWarning)
            settings_ase['outdir'] = self.workdir
            settings_ase['pseudo_dir'] = os.path.join(self.workdir, "pseudo")

        # Setting the new incut
        settings_ase['ecutwfc'] = self.encut

        # Saving updated pw.in
        with open(self.input_file, 'w') as f:
            espresso.write_espresso_in(f, atoms, input_data = settings_ase, kpts = (kpoint, kpoint, kpoint), pseudopotentials=pseudo)
            

    def calc(self) -> None:
        """Runs point simulation in the current folder
            I use my own compiled quantum espresso for that
        """
        output_file = os.path.join(self.workdir, f'{self.input_file_name.split(".")[0]}.out')
        self.calculator.run_scf(self.input_file, output_file)

    def extract_target(self, target:str) -> float:
        """Extracts the target from the simulation output"""
        out_file = os.path.join(self.workdir, f'{self.input_file_name.split(".")[0]}.out')
        data = read(out_file)
        if target == 'total_energy':
            return data.get_total_energy()