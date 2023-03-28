from enum import Enum


class PeriodicDftPackages(Enum):
    VASP = "vasp"
    vasp = "vasp"
    Vasp = "vasp"
    qe = "qe"
    QE = "qe"
    QuantumEspresso = "qe"
    quantumespresso = "qe"


class ConvergenceProperty(Enum):
    etotal = "total_energy"
    force = "force"
    phonon_modes = "phonon_modes"


class ConvergenceParameter(Enum):
    kpoints = "kpoints"
    encut = "encut"


def graceful_exit():
    print("aborting program")
    exit()


def messages(key):
    _msg = {
        "no_incar": "No {:s} exists at path {:s} specified",
        "multiple_incar": "Multiple matching {:s} files found, picking one of them",
        "no_incar_final": "No valid VASP incar files found",
        "no_calculation_input": "Either path to an INCAR or a _calculator object must be supplied",
        "set_to_singlepoint": "setting to a single point calculation for convergence tracking",
        "ignore_kspacing": "Ignoring KSPACING in input",
        "no_etotal": "Total energy not available yet, call VaspCalculation.run first",
        "multiple_structure_files": "more than one supported files found! Reading the first one: {:s}",
        "multiple_structure_images": "structure file {:s} has multiple structures, picking the first one",
        "converged": "Convergence found",
        "unconverged": "No convergennce found",
        "unconverged_kpts": "No converged kpoint mesh found!",
        "converged_kpts": "Converged k-point set: {:s}",
        "unsupported_package": "{:s} is not a supported periodic DFT package",
        "unsupported_property": "{:s} is not a supported convergence property",
        "unsupported_parameter": "{:s} is not a supported convergence parameter",
        "invalid_kpts": "invalid kpoints supplied: {:s}",
    }
    return _msg[key]
