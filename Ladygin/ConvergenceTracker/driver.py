from ase.io import read
import subprocess as sp

class Driver():
    """Base class for running simulation"""

    def __init__(self, input_file: str) -> None:
        """input_file - simulation settings"""
        self.input_file = input_file
        pass


    def gen_input(self, encut: float, kpoint: float):
        """Generates the input file for the driver based on input parameters
        
            encut - energy cutoff (in Ry)
            kpoint - kpoint dimension    
        """
        pass

    def run(self) -> None:
        """Runs point simulation in the current folder"""
        pass


    def extract_target(self, target:str) -> float:
        """Extracts the target from the simulation output"""


class DriverQE(Driver):

    def __init__(self, input_file: str = 'pw.in') -> None:
        super().__init__(input_file)

    def gen_input(self, encut: float, kpoint: float):

        with open(self.input_file, 'w') as f:
            f.write(f""" &control
    calculation     = 'scf'
    prefix          = 'si'
    restart_mode    = 'from_scratch'
    wf_collect      = .true.
    pseudo_dir      = './'
    outdir          = './'
    tprnfor         = .true.
    tstress         = .true.
 /
 &system
    ibrav           = 2
    celldm(1)       = 10.262 
    nat             = 2
    ntyp            = 1
    ecutwfc         = {encut}
    nbnd            = 8
 /
 &electrons
    diagonalization = 'david'
    mixing_beta     = 0.7
    conv_thr        = 1.0d-13
 /

ATOMIC_SPECIES
  Si  28.0855     Si-PBE.upf
ATOMIC_POSITIONS crystal
Si    0.00  0.00   0.00
Si   -0.25  0.75  -0.25
K_POINTS automatic
{kpoint} {kpoint} {kpoint} 0 0 0""")

    def run(self) -> None:

        sp.run("srun -n 4 -c 2 /pscratch/sd/v/vladygin/tools/q-e_new/q-e/bin/pw.x -nk 4 -input pw.in > pw.out", shell = True)

    def extract_target(self, target:str) -> float:
        data = read('pw.out')
        if target == 'total_energy':
            return data.get_total_energy()



if __name__ == '__main__':
    driver = DriverQE()

    driver.gen_input(40, 4)
    driver.run()
    total_energy = driver.extract_target('total_energy')
    print("Total energy is equal to", total_energy, 'eV')
