
import os
import numpy as np
import matplotlib.pyplot as plt
out_files = os.listdir('./results/')

total_energy = []
e_fermi = []
for prefix in out_files:
    temp_file = open(f'./results/{prefix}/scf.out', 'r').readlines()
    for i in range(len(temp_file)):
        if len(temp_file[i].split())>2:
            if temp_file[i].split()[0]=='!':
                temp_en=temp_file[i].split()[4]
            # if temp_file[i].split()[0]=='the':
            #     temp_fermi=temp_file[i].split()[-2]
    total_energy.append([float(temp_en),float(prefix)])
    # e_fermi.append([float(temp_fermi),float(file[:-13])])
    
tot_en = np.array(total_energy)
tot_en = tot_en[tot_en[:, 1].argsort()]


# e_f = np.array(e_fermi)
# e_f = e_f[e_f[:, 1].argsort()]

# fig = plt.figure(figsize=(8,6))
plt.plot(tot_en.T[1],tot_en.T[0])
# plt.ylabel('Total Energy [Ry]',fontsize=20)
# plt.xlabel('Smearing Value [Ry]',fontsize=20)
# plt.xticks(fontsize=20)
# plt.yticks(fontsize=20)
plt.savefig('total_sigma.png')
