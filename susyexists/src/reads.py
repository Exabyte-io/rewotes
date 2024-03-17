import numpy as np
import json

def read_pp(paths):
    elements = []
    wavefunction = []
    chargedensity = []
    for path in paths:
        pp = {"Element":"","Wavefunction":50,"Chargedensity:":200}
        with open(path, 'r') as data:
            data = data.read().split()
            for i in range(200):
                if data[i] =="Element:":
                    # print(str(f"Element: {data[i+1]}"))
                    pp["Element"]=data[i+1]
                if data[i] =="wavefunctions:":
                    # print(f"Wave function cutoff: {float(data[i+1])}")
                    pp["Wavefunction"]=data[i+1]
                if data[i] =="density:":
                    # print(f"Charge density cutoff: {float(data[i+1])}")
                    pp["Chargedensity"]=data[i+1]

        elements.append(pp['Element'])
        wavefunction.append(float(pp['Wavefunction']))
        chargedensity.append(float(pp['Chargedensity']))
    chargedensity.sort()
    wavefunction.sort()
    return(elements,wavefunction[-1]*1.5,chargedensity[-1]*1.5)


def read_poscar(path):
    data = open(path, 'r').readlines()
    cell =[]
    for i in data[2:5]:
        line=i.split()
        cell.append(line)
    atoms= []
    for i in data[8:]:
        line=(i.split())
        position=line[:3] 
        atom = [line[-1]]
        atoms.append(atom+position)
    return(cell,atoms)
    
        
def read_vc_relax(path):
    vc_relax = open(path, 'r').readlines()
    begin = 0
    end = 0
    for i in range(len(vc_relax)):
        if vc_relax[i] == 'Begin final coordinates\n':
            begin = i
        if vc_relax[i] == 'End final coordinates\n':
            end = i
    cell = np.array([i.split()
                    for i in vc_relax[begin+5:begin+8]]).astype(float)
    atoms = np.array([i.split() for i in vc_relax[begin+10:end]])
    return(cell, atoms)

def read_relax(path):
    data = open(path, 'r').readlines()
    begin = 0
    end = 0
    for i in range(len(data)):
        if data[i] == 'Begin final coordinates\n':
            begin = i
        if data[i] == 'End final coordinates\n':
            end = i
    atoms = np.array([i.split() for i in data[begin+3:end]])
    return atoms

def read_symmetries(fstring):
    f = open(fstring, 'r')
    x = np.zeros(0)
    for i in f:
        if "high-symmetry" in i:
            x = np.append(x, float(i.split()[-1]))
    f.close()
    return x

def read_efermi(path):
    lines = open(path, 'r').readlines()
    e_fermi = 0
    for i in lines:
        if "the Fermi energy is" in i:
            e_fermi = float(i.split()[-2])
            return e_fermi

def read_num_bands(self):
    directory = f"./Projects/{self.project_id}/{self.job_id}/scf.out"
    lines = open(directory, 'r').readlines()
    number = 0
    for i in lines:
        if "    number of Kohn-Sham states" in i:
            number = float(i.split()[-1])
            return number
        

def read_json(data=False,path=False):
    if (path):
        with open(path) as f:
            data = f.read()
    crystal = json.loads(data)
    atoms = []
    cell = np.array(crystal['lattice']['matrix'],dtype=str)
    for i in range(len(crystal['sites'])):
        element = crystal['sites'][i]['species'][0]['element']
        coordinate = crystal['sites'][i]['abc']
        atoms.append([element,str(coordinate[0]),str(coordinate[1]),str(coordinate[2])])
        # atoms.append(element)
        # coordinate = crystal['sites'][i]['abc']
        # positions[i]={element:coordinate}
    return cell,atoms
 
