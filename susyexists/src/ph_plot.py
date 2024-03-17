prefix = "0.0057"
labels = ['Γ',"M","K",'Γ']


import numpy as np
import matplotlib.pyplot as plt

def Symmetries(fstring):
    f = open(fstring, 'r')
    x = np.zeros(0)
    for i in f:
        if "high-symmetry" in i:
            x = np.append(x, float(i.split()[-1]))
    f.close()
    return x

freq = np.loadtxt("./{}.freq.gp".format(prefix))
load_sym = Symmetries("./plotband.out")
sym = load_sym/max(load_sym)
ph_path = freq.T[0]/max(freq.T[0])
cm2mev = 0.124

fig = plt.figure(figsize=(7,6))
# colors = ["blue","orange","green","red","purple","teal"]
for i in sym[:-1]:
    plt.axvline(i,linestyle="--",color="black")
for i in range(1,len(freq.T)):
    plt.plot(ph_path,freq.T[i]*cm2mev,linewidth=2,color="blue")
plt.xticks(sym,labels,fontsize=15)
# plt.axhline(0,color="black")
plt.title(f"{prefix}")
plt.ylim(0,)
plt.xlim(0,ph_path[-1])
plt.ylabel("ω (meV)",fontsize=15)
plt.savefig('phonon_band.png')
plt.show()