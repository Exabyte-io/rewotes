import argparse
import numpy as np
import matplotlib.pyplot as plt


ap = argparse.ArgumentParser()
ap.add_argument("-n", "--name", required=True,
                help="Project name")
ap.add_argument("-d", "--degauss", required=True, help="Degauss value")
args = vars(ap.parse_args())

project_name = args['name']
degauss = args['degauss']


def read_efermi(path):
    lines = open(path, 'r').readlines()
    e_fermi = 0
    for i in lines:
        if "the Fermi energy is" in i:
            e_fermi = float(i.split()[-2])
            return e_fermi
        
def Symmetries(fstring):
    f = open(fstring, 'r')
    x = np.zeros(0)
    for i in f:
        if "high-symmetry" in i:
            x = np.append(x, float(i.split()[-1]))
    f.close()
    return x

sym = Symmetries(f'./{project_name}/{degauss}/bands-pp.out')


temp_file = open(f'./{project_name}/{degauss}/scf.out', 'r').readlines()
for i in range(len(temp_file)):
    if len(temp_file[i].split())>2:
        if temp_file[i].split()[0]=='the':
            fermi=temp_file[i].split()[-2]

fig = plt.figure(figsize=(8,6))
data =np.loadtxt(f'./{project_name}/{degauss}/bands.dat.gnu')

k = np.unique(data[:, 0])
bands = np.reshape(data[:, 1], (-1, len(k)))
for band in range(len(bands)):
    plt.plot(k, bands[band, :]-float(fermi),c='black')
plt.xticks(sym,['G','M','K','G'])
plt.axvline(sym[1],c='black')
plt.axvline(sym[2],c='black')
plt.axhline(0,c='red')
plt.ylim(-10,10)
plt.xlim(sym[0],sym[-1])
plt.savefig(f'./{project_name}/plots/band_{degauss}.png')
# plt.show()
                