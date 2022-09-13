from computing.executable import Executable
from crystal import reciprocal
from crystal.voronoi import Voronoi
from simulation.qe_spec import QESpec
from simulation.qe_process import QEProcess
from simulation.simulation import Simulation

import math
import os
import numpy as np

class Convergence(Executable):
    """
    Class responsible for running a series of
    Simulations to do a convergence test for
    a given material system.
    """
    
    def __init__(self, dir, encut, thresh, rsx, input_name="input.txt", meta=None, obtusify=True, homogeneous_k=False):
        """
        Sets up a convergence test in the specified
        directory, using the provided QESpec, and
        a metadata string for logging purposes.
        """
        self.dir = dir
        self.spec = QESpec.parse_file(os.path.join(dir, input_name))
        self.spec.set_encut(encut)
        self.thresh = thresh
        if obtusify:
            Convergence._set_obtuse_reciprocal_lattice(self.spec)
        lattice = self.spec.get_lattice()
        r_lattice = np.transpose(reciprocal.get_reciprocal(lattice))
        r_lengths = r_lattice[0]*r_lattice[0]
        for i in range(1, len(r_lattice)):
            r_lengths += r_lattice[i]*r_lattice[i]
        self.r_lengths = list(r_lengths)
        for i in range(0, len(self.r_lengths)):
            self.r_lengths[i] = self.r_lengths[i]**0.5
        self.meta = meta
        self.homogeneous_k = homogeneous_k
        self.rsx = rsx
        self.k = self.spec.get_k_points()
        self.ks = []
        self.values = []
      
    @staticmethod
    def _set_obtuse_reciprocal_lattice(spec):
        if spec.get_lattice_type() != 0:
            return
        try:
            start = spec.get_lattice()
            reci = reciprocal.get_reciprocal(start)
            reci_obtuse = Voronoi.obtuse_basis(reci)
            real_reci_obtuse = reciprocal.get_reciprocal(reci_obtuse)
            unimodular = np.matmul(np.linalg.inv(start), real_reci_obtuse)
            b_start = spec.get_positions()
            b_end = list(np.matmul(np.linalg.inv(unimodular), b_start[1]))
            spec.set_lattice(real_reci_obtuse)
            spec.set_positions(b_start[0], b_end)
        except Exception, e:
            return
    
    def get_results(self, analysis):
        """
        Extract the value used to test convergence using the
        provided QEAnalysis object.
        """
        pass
    
    def get_resources(self):
        return self.rsx
        
    def run(self, envr):
        sim_start = self.next_executable()
        if not sim_start.run(envr):
            print("Error encountered in starting simulation.")
            return False
        self.next_k()
        sim_next = self.next_executable()
        if not sim_next.run(envr):
            print("Error encountered in simulation.")
            return False
        while abs(self.values[len(self.values) - 1] -
                  self.values[len(self.values) - 2]) > self.thresh:
            self.next_k()
            sim_next = self.next_executable()
            if not sim_next.run(envr):
                print("Error encountered in simulation.")
                return False
        
    def next_executable(self):
        self.spec.set_k_points(self.k)
        process = Convergence.CustomProcess(self.dir, self.spec, self)
        sim = Simulation.construct_override(self.dir, self.spec, self.rsx, None, None, process)
        return sim
    
    def next_k(self):
        last_k = self.k
        if self.homogeneous_k:
            self.k = [x + 1 for x in last_k]
        else:
            lin_density = [last_k[i]/self.r_lengths[i] for i in range(0, len(self.r_lengths))]
            choice = lin_density.index(min(lin_density))
            self.k = [x for x in last_k]
            self.k[choice] += 1
        
    class CustomProcess(QEProcess):
        def __init__(self, dir, spec, parent):
            super(Convergence.CustomProcess, self).__init__(dir, spec)
            self.parent = parent
        def get_results(self, analysis):
            outcome = self.parent.get_results(analysis)
            self.parent.ks.append(self.parent.k)
            self.parent.values.append(outcome)
            return outcome
