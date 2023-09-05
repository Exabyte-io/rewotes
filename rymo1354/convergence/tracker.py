import numpy as np
from copy import deepcopy
from tabulate import tabulate
from ase.optimize import QuasiNewton
import sys

class GenerateKpoints(object):
    def __init__(self):
        pass

    def my_ceil(self, value, precision=6):
        ''' Deals with floating point uncertainty '''
        return np.true_divide(np.ceil(value * 10**precision), 10**precision)

    def calc_Kspacing(self, reciprocal_vector, kpts):
        return np.linalg.norm(reciprocal_vector)*(2*np.pi)/kpts

    def calc_N(self, reciprocal_vector, kspacing):
        return np.max([1, int(np.ceil(np.linalg.norm(reciprocal_vector)*((2*np.pi)/kspacing)))])

    def calc_new_division(self, reciprocal_lattice, old_division):
        increase_kspacings = [self.my_ceil(self.calc_Kspacing(reciprocal_lattice[0], old_division[0] + 1)), 
                              self.my_ceil(self.calc_Kspacing(reciprocal_lattice[1], old_division[1] + 1)), 
                              self.my_ceil(self.calc_Kspacing(reciprocal_lattice[2], old_division[2] + 1))]
        new_kspacing = np.max(increase_kspacings) # Largest kspacing that induces a change in kpts
        new_division = [self.calc_N(reciprocal_lattice[0], new_kspacing), 
                        self.calc_N(reciprocal_lattice[1], new_kspacing),
                        self.calc_N(reciprocal_lattice[2], new_kspacing)] 
        return new_division

class ConvergenceHandler(object):
    def __init__(self, atoms, calculator, kpoints, threshold, condition='subsequent', prop='total energy'):
        self.atoms = atoms
        self.calculator = calculator
        if 'ediffg' in self.calculator.asdict()['inputs']:
            if np.sign(self.calculator.asdict()['inputs']['ediffg']) == -1:
                self.fmax = np.absolute(self.calculator.asdict()['inputs']['ediffg'])
            else:
                self.fmax = 1000 # Only total energy converged 
        else:
            self.fmax = 1000 # Only total energy converged
        self.starting_kpoints = kpoints
        self.threshold = threshold
        self.condition = condition
        self.prop = prop
        self.generate_kpoints = GenerateKpoints()

    def convergence_property(self, opt_atoms):
        if self.prop == 'total energy':
             return opt_atoms.get_potential_energy() 
        else: # to support other properties
             print('No other properties supported at this time')
             sys.exit(1)
    
    def convergence_condition(self, all_props):
        if self.condition == 'subsequent':
            if len(all_props) >= 2:
                old_prop, new_prop = all_props[-2], all_props[-1]
                abs_diff = np.absolute(np.subtract(old_prop, new_prop))
                if abs_diff <= self.threshold:
                    return True
                else:
                    return False
            else:
                return False
        else:
             print('No other convergence conditions supported at this time')
             sys.exit(1)

    def make_table(self, all_kpoints, all_props, all_diffs):
        table = []
        for i in range(len(all_kpoints)):
            table.append([all_kpoints[i], all_props[i], all_diffs[i]])
        return tabulate(table, headers=['KPOINTS', self.prop, 'difference'])

    def converge_kpoints(self):
        all_kpoints = []
        all_props = []
        all_diffs = []
        kpoints = deepcopy(self.starting_kpoints)
        
        while self.convergence_condition(all_props) == False:
            self.calculator.set(kpts=kpoints)
            opt_atoms = deepcopy(self.atoms)
            opt_atoms.set_calculator(self.calculator)

            kspacing = np.min(self.generate_kpoints.calc_Kspacing(self.atoms.get_reciprocal_cell(), kpoints))
            r_kspacing = self.generate_kpoints.my_ceil(kspacing)

            print('Running with %s KPOINTS using KSPACING = %s Angstroms-1' % (str(kpoints), str(r_kspacing)))
            try:
                optimizer = QuasiNewton(opt_atoms)
                optimizer.run(fmax=self.fmax)
                all_kpoints.append(kpoints)
                all_props.append(opt_atoms.get_potential_energy())
                try:
                    all_diffs.append(np.absolute(all_props[-2]-all_props[-1]))
                except IndexError:
                    all_diffs.append('N/A')
            except IndexError:
                # Can occur if too few KPOINTS for tetrahedron method
                print('Bad calculation, check input files')
            kpoints = self.generate_kpoints.calc_new_division(self.atoms.get_reciprocal_cell(), kpoints) 
 
        final_kpoints = str(all_kpoints[-2])
        convergence_table = self.make_table(all_kpoints, all_props, all_diffs)
        
        return final_kpoints, convergence_table 
    


