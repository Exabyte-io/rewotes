import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import re
import subprocess
import json
import qcelemental as qcel
import untangle


def make_monolayer(atoms):
    df = pd.DataFrame()
    atoms=np.array(atoms)
    df['Atoms'] = atoms.T[0]
    df['x'] = atoms.T[1].astype(float)
    df['y'] = atoms.T[2].astype(float)
    df['z'] = atoms.T[3].astype(float)
    df = df.query('z<0.5')
    df_shifted = pd.DataFrame()
    df_shifted["Atoms"] = df['Atoms'].values
    df_shifted["x"] = df['x'].values
    df_shifted["y"] = df['y'].values
    df_shifted["z"] = df['z'].values+0.25
    return df_shifted.values


def plot_sigma_energy(path):
    out_files = os.listdir(path)
    total_energy = []
    e_fermi = []
    for file in out_files:
        temp_file = open(f'{path}/{file}', 'r').readlines()
        for i in range(len(temp_file)):
            if len(temp_file[i].split())>2:
                if temp_file[i].split()[0]=='!':
                    temp_en=temp_file[i].split()[4]
                # if temp_file[i].split()[0]=='the':
                #     temp_fermi=temp_file[i].split()[-2]
        total_energy.append([float(temp_en),float(file[:-13])])
        # e_fermi.append([float(temp_fermi),float(file[:-13])])
        
    tot_en = np.array(total_energy)
    tot_en = tot_en[tot_en[:, 1].argsort()]
    
    
    # e_f = np.array(e_fermi)
    # e_f = e_f[e_f[:, 1].argsort()]
    
    fig = plt.figure(figsize=(8,6))
    plt.plot(tot_en.T[1],tot_en.T[0])
    plt.ylabel('Total Energy [Ry]',fontsize=20)
    plt.xlabel('Smearing Value [Ry]',fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.savefig('total_sigma.png')

def get_total_energy(path):
    obj = untangle.parse(path)
    en = float(obj.qes_espresso.output.total_energy.etot.cdata)*2
    # p = subprocess.Popen(f"grep '!' {path}", shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    # line = p.stdout.readlines()[0].decode()
    # en = float(re.findall(r"[-+]?(?:\d*\.*\d+)", line)[0])
    return(en)

def configure(calculation,path="./config.json"):
    with open(path) as f:
        data = f.read()
        config = json.loads(data)
    return config[calculation]

def afm_maker(atom,afm_matrix):
    # afm_matrix = [['u','u','d','d'],['u','d','u','d'],['u','d','d','u']]
    afm =[]
    for i in range(len(afm_matrix)):
        # print(f"AFM{i}")
        afm.append(np.array(atom).copy())
        for j in range(4):
                afm[i][j][0]=atom[j][0]+afm_matrix[i][j]
                # print(afm[i][j][0])
    return afm

def default_pseudo(atom):
    atom_type = list(set([a[0] for a in atom]))
    atom_array = []
    for i in atom_type:
        temp_atom = {"atom":i,'mass':str(qcel.periodictable.to_mass(i)),'pseudopotential':f"{i}.UPF"}
        atom_array.append(temp_atom)
        # print(temp_atom)
    # print(atom_array)  
    return atom_array  

def atom_type(atom):
    num_type = len(list(set([a[0] for a in atom])))
    return num_type


def test_parameter(self,parameter_name,conv_thr,start,end,step,num_core,debug=False,out=False):
    parameter = np.arange(start,end,step)
    result = np.zeros(shape=(3,len(parameter)))
    end=0
    for j,i in enumerate(parameter):
        if parameter_name=="ecutwfc":
            self.ecutwfc(i)
            self.job_id=f"ecutwfc_{i}"
        if parameter_name=="kpoints":
            self.k_points(int(i))
            self.job_id=f"kpoints_{i}"
        if debug==False:
            self.scf(num_core)
        path = f'./Projects/{self.project_id}/{self.job_id}/{self.job_id}.save/data-file-schema.xml'
        temp_en = get_total_energy(path)
        temp_time = get_time(path)
        result[0][j]=i #parameters
        result[1][j]=temp_en #energy
        result[2][j]=temp_time #time
        end+=1
        if j!=0:
            if out==True:
                print(f"{parameter_name}: {result[0][j]}     DeltaE :{(result[1][j]-result[1][j-1])} Ry    Time: {result[2][j]} seconds")
            if abs(result[1][j-1] - result[1][j]) < conv_thr:
                end=j
                break
    # print(result.T[:end+1].T)
    return result.T[:end+1].T


def get_time(path):
    obj = untangle.parse(path)
    time = float(obj.qes_espresso.timing_info.total.wall.cdata)
    return time
    # with open(path, 'r') as data:
    #  data = data.read().split()
    #  counter = 0 
    #  for j,i in enumerate(data):
    #     # if i == "PWSCF":
    #     #     counter += 1 
    #     #     if counter ==3:
    #     #                         # print(f"CPU: {data[j+2]}, WALL: {data[j+4]}")
    #     #                         return (data[j+4])
    #     if i=='End':
    #        return(data[j-2])
