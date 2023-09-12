import numpy as np
from copy import deepcopy
from tabulate import tabulate
from ase.optimize import QuasiNewton
from convergence.constructors import CalculatorConstructor, AtomsConstructor
from convergence.constructors import GetFmax, PropertyConstructor
import sys
import os

class KptsOperator:
    def __init__(self):
        pass

    def my_ceil(self, value, precision=6):
        ''' Deals with floating point uncertainty '''
        return np.true_divide(np.ceil(value * 10**precision), 10**precision)
 
    def calc_Kspacing(self, reciprocal_vector, kpts):
        return np.linalg.norm(reciprocal_vector)*(2*np.pi)/kpts

    def calc_N(self, reciprocal_vector, kspacing):
        return np.max([1, int(np.ceil(np.linalg.norm(reciprocal_vector)*((2*np.pi)/kspacing)))])

    def calc_new_division(self, atoms, old_division):
        reciprocal_lattice = atoms.get_reciprocal_cell()
        new_kspacing = np.max([self.my_ceil(self.calc_Kspacing(reciprocal_lattice[i], 
                               old_division[i] + 1)) for i in range(len(reciprocal_lattice))])
        new_division = [self.calc_N(reciprocal_lattice[i], new_kspacing) for i in range(len(reciprocal_lattice))]
        return new_division

    def calc_used_spacing(self, atoms, division):
        reciprocal_lattice = atoms.get_reciprocal_cell()
        used_spacing = np.max([self.my_ceil(self.calc_Kspacing(reciprocal_lattice[i],
                               division[i])) for i in range(len(reciprocal_lattice))])
        return used_spacing

class KptsConvergenceTracker():
    def __init__(self, convergence_property, convergence_threshold):
        self.all_kpts = []
        self.properties = []
        self.differences = []
        
        self.convergence_threshold = convergence_threshold
        self.convergence_property = convergence_property

    def __converged(self):
        try:
            if self.differences[-1] <= self.convergence_threshold:
                return True
            else:
                return False
        except (IndexError, TypeError): #No numerical difference here
            return False

    def __initialize_run(self, calculation_type, directory, **kwargs):
        atoms = AtomsConstructor().make(key=calculation_type, directory=directory)
        calculator = CalculatorConstructor().make(key=calculation_type, directory=directory, **kwargs)
        fmax = GetFmax().make(key=calculation_type, calculator=calculator)
        get_Kpts = KptsOperator()
        return atoms, calculator, fmax, get_Kpts

    def __converge_kpts(self, atoms, calculator, get_Kpts, fmax, kpoints):
        while self.__converged() == False:
            print('Running %s, KSPACING = %s A-1' % (str(kpoints), str(get_Kpts.calc_used_spacing(atoms, kpoints))))
            calculator.set(kpts=kpoints)
            opt_atoms = deepcopy(atoms)
            opt_atoms.set_calculator(calculator)

            try:
                optimizer = QuasiNewton(opt_atoms)
                optimizer.run(fmax=fmax)
                try:
                    self.all_kpts.append(kpoints)
                    self.properties.append(PropertyConstructor().make(key=self.convergence_property, atoms=opt_atoms))
                    self.differences.append(np.absolute(np.subtract(self.properties[-2], self.properties[-1])))
                except IndexError:
                    self.differences.append('N/A')
            except IndexError:
                # Can occur if too few KPOINTS for tetrahedron method
                print('Bad calculation, check input files')
            kpoints = get_Kpts.calc_new_division(atoms, kpoints)
        return

    def __make_table(self):
        table = []
        for i in range(len(self.all_kpts)):
            table.append([self.all_kpts[i], self.properties[i], self.differences[i]])
        return tabulate(table, headers=['KPOINTS', self.convergence_property, 'difference'])

    def run(self, calculation_type, directory, **kwargs):
        atoms, calculator, fmax, get_Kpts = self.__initialize_run(calculation_type, directory, **kwargs)
        self.__converge_kpts(atoms, calculator, get_Kpts, fmax, kpoints=[1, 1, 1])
        final_kpoints = str(self.all_kpts[-2])
        convergence_table = self.__make_table()
        return final_kpoints, convergence_table
