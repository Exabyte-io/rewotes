from ase import Atoms
import os
# Some usefull consts
consts = {'Ry': 13.6}

# Silicon template in case no pw.in
pw_template = {'control': {'calculation': 'scf',
  'prefix': 'si',
  'restart_mode': 'from_scratch',
  'wf_collect': '.true.',
  'pseudo_dir': './pseudo',
  'outdir': './',
  'tprnfor': '.true.',
  'tstress': '.true.'},
 'system': {'ibrav': 2.0,
  'celldm(1)': 10.262,
  'nat': 2.0,
  'ntyp': 1.0,
  'ecutwfc': 40.0,
  'nbnd': 8.0},
 'electrons': {'diagonalization': 'david',
  'mixing_beta': 0.7,
  'conv_thr': '1.0d-13'},
 'ATOMIC_SPECIES': {'Si': '28.0855 Si-PBE.upf'}}

atoms_template = Atoms(symbols='Si2', pbc=True, cell=[[-2.7152082560620205, 0.0, 2.7152082560620205], [0.0, 2.7152082560620205, 2.7152082560620205], [-2.7152082560620205, 2.7152082560620205, 0.0]], positions = [[0.        , 0.        , 0.        ], [1.35760413, 1.35760413, 1.35760413]])


def parse_qe_input(path: str) -> dict:
    """
    Parse the input file of Quantum Espresso and return it as a dictionary
    of the parameters.

    Args:
        path (:obj:`str`):
            Path to the Quantum Espresso input file.

    Returns: :obj:`dict`
    """
    if not os.path.exists(path):
        raise FileNotFoundError(f"File {path} does not exist!")
    sections = [
        "ATOMIC_SPECIES",
        "ATOMIC_POSITIONS",
        "K_POINTS",
        "ADDITIONAL_K_POINTS",
        "CELL_PARAMETERS",
        "CONSTRAINTS",
        "OCCUPATIONS",
        "ATOMIC_VELOCITIES",
        "HUBBARD",
        "ATOMIC_FORCES",
        "SOLVENTS",
    ]
    with open(path, encoding="utf-8") as f:
        lines = f.readlines()
    params = {}
    for raw_line in lines:
        line = raw_line.strip().partition("#")[0]  # ignore in-line comments
        if line.startswith("&") or any(sec in line for sec in sections):
            section = line.split()[0].replace("&", "")
            params[section] = {}
        elif line.startswith(("#", "!", "/")):
            continue
        elif "=" in line:
            key, value = line.split("=")
            # Convent numeric values to float
            try:
                value = float(value.strip())
            except Exception:
                # string keywords in QE input file are enclosed in quotes
                # so we remove them to avoid too many quotes when generating
                # the input file with ase
                value = value.strip()
                value = value.replace("'", "")
                value = value.replace('"', "")

            params[section][key.strip()] = value
        elif len(line.split()) > 1:
            key, value = line.split()[0], " ".join([str(val) for val in line.split()[1:]])
            params[section][key.strip()] = value
    # Remove structure info (if present), as will be re-written with distorted structures
    for section in [
        "ATOMIC_POSITIONS",
        "K_POINTS",
        "ADDITIONAL_K_POINTS",
        "CELL_PARAMETERS",
    ]:
        params.pop(section, None)
    if "SYSTEM" in params:
        for key in ["celldm(1)", "nat", "ntyp", "ibrav"]:
            params["SYSTEM"].pop(key, None)
    return params


def convert_settings(settings:dict)->[dict, dict]:
    """Converts qe settings to ase format
    
        returns settings and pseudo in ase format 
    """
    new_dict = {key: settings[tag][key] for tag in ['CONTROL', 'SYSTEM', 'ELECTRONS', 'control', 'system', 'electrons'] if tag in settings for key in settings[tag]}

    for key in new_dict:
        # Fixing floats that should be ints
        try:
            if int(new_dict[key]) == new_dict[key]:
                new_dict[key] = int(new_dict[key])
        except:
            pass
        # Fixing exponential format
        try:
            new_dict[key] = float(new_dict[key].replace('d', 'e'))
        except:
            pass
    
    for key in new_dict:
        if new_dict[key] == '.true.':
            new_dict[key] = True
        if new_dict[key] == '.false.':
            new_dict[key] = False
    
    pseudo = {key: settings["ATOMIC_SPECIES"][key].split()[1] for key in settings["ATOMIC_SPECIES"]}
    return new_dict, pseudo