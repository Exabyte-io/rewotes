import numpy as np
import matplotlib.pyplot as plt
from . import reads

def plot_electron(self,ylim=False,show=False,save=False):
    sym = reads.read_symmetries(f'./Projects/{self.project_id}/{self.job_id}/bands-pp.out')
    fermi = reads.read_efermi(f'./Projects/{self.project_id}/{self.job_id}/scf.out')
    fig = plt.figure(figsize=(8,6))
    data = np.loadtxt(f'./Projects/{self.project_id}/{self.job_id}/bands.dat.gnu')
    k = np.unique(data[:, 0])
    bands = np.reshape(data[:, 1], (-1, len(k)))
    for band in range(len(bands)):
        plt.plot(k, bands[band, :]-float(fermi),c='black')
    plt.xticks(sym,self.path)
    for i in range(1,len(sym)-1):
        plt.axvline(sym[i],c='black')
    plt.axhline(0,c='red')
    plt.ylim(ylim[0],ylim[1])
    plt.xlim(sym[0],sym[-1])
    if save==True:
        plt.savefig(f'./Projects/{self.project_id}/{self.job_id}/band.png')
    if show==True:
        return fig
    

def plot_phonon(self):
    sym = []
    point = [0]
    for k,i in enumerate(self.config['k_points_bands']):
        sym.append(i['label'].split()[1])
        if k!=len(self.config['k_points_bands'])-1:
            point.append(point[k]+int(i['number']))
    freq = np.loadtxt(f"./Projects/{self.project_id}/{self.job_id}/{self.job_id}.freq.gp")
    ph_path = freq.T[0]/max(freq.T[0])
    cm2mev = 0.124
    fig = plt.figure(figsize=(7,6))
    
    for i in range(1,len(freq.T)):
        plt.plot(ph_path,freq.T[i]*cm2mev,linewidth=2,color="blue")
    # print(point)
    tick = [ ph_path[i] for i in point ]
    for i in tick[1:-1]:
        plt.axvline(i,linestyle="--",color="black")
    plt.xticks(tick,sym,fontsize=15)
    plt.ylim(0,)
    plt.xlim(0,ph_path[-1])
    plt.ylabel("ω (meV)",fontsize=15)
    plt.savefig('phonon_band.png')
    plt.show()
    
