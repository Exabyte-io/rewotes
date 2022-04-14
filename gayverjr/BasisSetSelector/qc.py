import os
import tempfile
import subprocess
from io import StringIO
import traceback

# https://physics.nist.gov/cgi-bin/cuu/Value?hrev
AU2EV= 27.211386245988

# TODO: custom basis sets (specified by path to formatted file)
# TODO: other properties (e.q. vibrational frequenices)
# TODO: handling when NWChem is not installed

def get_nwchem_homo_lumo_gap(output):
    output_stream = StringIO(output)
    lines = output_stream.readlines()
    idx = 0
    while 'DFT Final Molecular Orbital Analysis' not in lines[idx]:
        idx+=1
    prev_E = None
    while idx<len(lines):
        if 'Occ=' in lines[idx]:
            split_line = lines[idx].replace('=','= ').split()
            occ = float(split_line[3].replace("D", "E"))
            if occ == 0.0:
                lumo_E = float(split_line[5].replace("D", "E"))
                return {'success':True,'result':[(lumo_E-prev_E)*AU2EV]}
            else:
                prev_E = float(split_line[5].replace("D", "E"))
        idx+=1
    return {'success':False}

def run_nwchem(basis,mol,functional,prop_type):
    work_dir = tempfile.TemporaryDirectory()
    with open(os.path.join(work_dir.name,'nwchem.nw'),'w') as f:
        tasks = ['task','dft']
        f.write('start Job{}\n'.format(mol.id))
        f.write('title \" Job{}\"\n'.format(mol.id))
        f.write('geometry units ang\n')
        f.write(str(mol)+'\n')
        f.write('end\n\n')
        f.write('basis\n')
        f.write('* library {}\n'.format(basis.name))
        f.write('end\n\n')
        f.write('dft\n')
        f.write('XC  {}\n'.format(functional))
        f.write('end\n\n')
        if prop_type=='homo lumo gap':
            pass
        for task in tasks:
            f.write('{} '.format(task))
    try:
        output = subprocess.run(['nwchem'], shell=True, capture_output=True,encoding='UTF-8',cwd=work_dir.name)
        work_dir.cleanup()
        return get_nwchem_homo_lumo_gap(output.stdout)
    except Exception as e:
        work_dir.cleanup()
        return {'success':False}
    

def run_qc(basis,mol,functional,prop_type,engine):
    if engine=='nwchem':
        return run_nwchem(basis,mol,functional,prop_type)
    else:
        raise RuntimeError("Only NWChem currently supported.")