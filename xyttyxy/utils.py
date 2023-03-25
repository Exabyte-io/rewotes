from enum import Enum


class PeriodicDftPackages(Enum):
    VASP = 'vasp'
    vasp = 'vasp'
    Vasp = 'vasp'
    qe = 'qe'
    QE = 'qe'
    QuantumEspresso = 'qe'
    quantumespresso = 'qe'
    

class ConvergenceProperty(Enum):
    etotal = 'total_energy'
    force = 'force'
    phonon_modes = 'phonon_modes'


class ConvergenceParameter(Enum):
    kpoints = 'kpoints'
    encut = 'encut'


def graceful_exit():
    print('aborting program')
    exit()
