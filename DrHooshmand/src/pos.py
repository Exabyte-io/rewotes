import numpy as np
import math
import ast
import sys
np.set_printoptions(16)

class POSCAR:
    """
    Reading POSCAR
    """

    def __init__(self, name = None, label = None, a = None, pos = None, pos_inv = None,
                chemistry = None, chem_lab = None, N = None, style = None, atoms = None):
        """
        Create an instance of POSCAR
        """
        if name is not None:
            self.selective, self.label, self.a, self.pos, self.chemistry, self.chem_lab, self.N, self.style, self.atoms = self.rPOSCAR(name)
            self.pos_inv = np.linalg.inv(self.pos)


    @staticmethod
    def rPOSCAR (name,thres=1e-6):

        selective = False
        f = open(name, 'r')
        label = f.readline()
        label = label.strip()
        a = f.readline()
        a = float(a.strip())

        pos = []
        for i in range(3):
            line = f.readline()
            line = line.strip()
            y = [float(j) for j in line.split()]
            pos.append(y)

        pos = np.array(pos)
        pos[abs(pos)<thres]=0.0
        pos.tolist()

        chem_lab = f.readline()
        chem_lab = chem_lab.strip().split()
        line = f.readline()
        line = line.strip()
        chemistry = [int(j) for j in line.split()]

        N = 0
        for i in chemistry:
            N += i

        y = f.readline()
        y = y.strip()
        if y[0]== "S":
            selective =True
            print("Selective Dynamics Read")
            style = f.readline()
            style = style.strip()
        else:
            style = y.strip()

        atoms = []
        for i in range(N):
            line = f.readline()
            line = line.strip()
            # print(line)
            if selective:
                y = [float(j) for j in line.split()[0:3]]
            else:
                y = [float(j) for j in line.split()]
            # print(y)
            atoms.append(y)

        atoms = np.array(atoms)
        atoms[abs(atoms)<thres] =0.0
        atoms.tolist()

        f.close()

        return selective, label, a, np.array(pos), chemistry, chem_lab, N, style, np.array(atoms)



    def wPOSCAR(self, name,label):

        f = open(name,'w')

        f.write(label + '\n')
        f.write(str(self.a) + '\n')

        for i in range(3):
            f.write(str(self.pos[i][0]) + " " + str(self.pos[i][1]) + " " + str(self.pos[i][2])+  '\n')

        for i in self.chem_lab:
            f.write(str(i)+" ")
        f.write('\n')

        for i in self.chemistry:
            f.write(repr(i)+" ")
        f.write('\n')
        f.write(self.style + '\n')

        for i in self.atoms:
            f.write(repr(i[0]) + " " + repr(i[1]) + " " + repr(i[2]) + '\n')
        f.close()
