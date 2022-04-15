import os
import tempfile
import subprocess
from io import StringIO
import traceback

# https://physics.nist.gov/cgi-bin/cuu/Value?hrev
AU2EV= 27.211386245988

# TODO: handling when NWChem is not installed
# TODO: support for mixed basis?
# TODO: logging?

def get_nwchem_frequencies(natom,output:str)->dict:
    ''' Parses frequencies from NWChem output.
    '''
    # Nwchem prints trivial modes, need to discard them
    output_stream = StringIO(output)
    lines = output_stream.readlines()
    idx = 0
    while idx<len(lines) and 'Normal Eigenvalue' not in lines[idx]:
        idx+=1
    if idx>=len(lines):
        return {'success':False}
    idx+=3
    nmodes = 3*natom
    freqs = []
    for i in range(0,nmodes):
        freqs.append(float(lines[idx].split()[1]))
        idx+=1
    freq_start = 6
    # NWChem kindly prints whether it is a linear molecule
    while idx<len(lines):
        if 'Linear Molecule' in lines[idx]:
            freq_start = 5
            break
        idx+=1
    return {'success':True,'result':freqs[freq_start:]}
    
def get_nwchem_homo_lumo_gap(output:str)->dict:
    ''' Parses homo lumo gap from NWChem output.
    '''
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
    ''' Constructs NWChem input and executes job as a subprocess.
    '''
    work_dir = tempfile.TemporaryDirectory()
    with open(os.path.join(work_dir.name,'nwchem.nw'),'w') as f:
        tasks = ['task','dft']
        f.write('start Job{}\n'.format(mol.id))
        f.write('title \" Job{}\"\n'.format(mol.id))
        f.write('geometry units ang\n')
        f.write(str(mol)+'\n')
        f.write('end\n\n')
        if any(pattern in basis.name for pattern in ['6-31','3-21','4-21','4-31']):
            f.write('basis cartesian\n')
        else:
            f.write('basis\n')
        f.write('* library {}\n'.format(basis.name))
        f.write('end\n\n')
        f.write('dft\n')
        f.write('XC  {}\n'.format(functional))
        f.write('end\n\n')
        if prop_type=='homo lumo gap':
            pass
        elif prop_type=="frequencies":
            tasks.append('frequencies')
            f.write('freq\nend\n\n')
        for task in tasks:
            f.write('{} '.format(task))
    try:
        output = subprocess.run(['nwchem'], shell=True, capture_output=True,encoding='UTF-8',cwd=work_dir.name)
        work_dir.cleanup()
        if prop_type == 'homo lumo gap':
            return get_nwchem_homo_lumo_gap(output.stdout)
        elif prop_type == 'frequencies':
            return get_nwchem_frequencies(mol.natoms,output.stdout)
        else:
            raise NotImplementedError("Only homo lumo gap and frequencies are implemented.")
    except Exception as e:
        work_dir.cleanup()
        return {'success':False}
    

def run_qc(basis,mol,functional,prop_type,engine):
    ''' Executes quantum chemistry calculation.
    '''
    if engine=='nwchem':
        return run_nwchem(basis,mol,functional,prop_type)
    else:
        raise RuntimeError("Only NWChem currently supported.")