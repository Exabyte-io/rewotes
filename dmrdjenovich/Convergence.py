from computing import Executable
from crystal import reciprocal
from crystal import voronoi

from abc import ABC, abstractmethod
import math
import os
import numpy as np

class Convergence(Executable):
    """
    Class responsible for running a series of
    Simulations to do a convergence test for
    a given material system.
    """
    
    INPUT_NAME = "input.txt"
    
    def __init__(self, dir, encut, thresh, rsx, meta=None, obtusify=True, homogeneous_k=False):
        """
        Sets up a convergence test in the specified
        directory, using the provided QESpec, and
        a metadata string for logging purposes.
        """
        self.dir = dir
        self.spec = QESpec.parse_file(os.join("/Users/david/Documents/Internship/Mat3ra/rewotes/dmrdjenovich/simulation", Convergence.INPUT_NAME))
        self.spec.set_encut(encut)
        self.thresh = thresh
        if obtusify:
            Convergence._set_obtuse_reciprocal_lattice(spec)
        lattice = spec.get_lattice()
        r_lattice = np.transpose(reciprocal.get_reciprocal(self.lattice))
        self.r_lengths = r_lattice[0]*r_lattice[0]
        for i in range(1, len(r_lattice)):
            self.r_lengths += self.r_lattice[i]*self.r_lattice[i]
        self.r_lengths = list(self.r_lengths)
        for i in range(0, len(self.r_lengths)):
            self.r_lengths[i] = self.r_lengths[i]**0.5
        self.meta = meta
        self.rsx = rsx
        self.k = spec.get_k_points()
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
    
    @abstractmethod
    def get_results(self, analysis):
        pass
    
    def get_resources(self):
        return self.rsx
        
    def exec(self, envr):
        sim_start = self.next_executable()
        if not sim_start.exec(envr):
            print("Error encountered in starting simulation.")
            return False
        self.next_k()
        sim_next = self.next_executable()
        if not sim_next.exec(envr):
            print("Error encountered in simulation.")
            return False
        while math.abs(self.vals[len(self.vals) - 1] - self.vals[len(self.vals) - 2]) > self.thresh:
            self.next_k()
            sim_next = self.next_executable()
            if not sim_next.exec(envr):
                print("Error encountered in simulation.")
                return False
        
    def next_executable(self):
        self.spec.set_k_points(self.k)
        process = Convergence.CustomProcess(self.dir, self.spec, self)
        sim = Simulation.construct_override(self.dir, self.spec, self.rsx, None, None, process)
        return sim
    
    def next_k(self):
        last_k = self.ks[len(self.ks) - 1]
        if self.homogeneous_k:
            self.k = [x + 1 for x in last_k]
        else:
            lin_density = [last_k[i]/self.r_lengths[i] for i in range(0, len(last_k))]
            choice = lin_density.index(min(lin_density))
            self.k = [x for x in last_k]
            self.k[choice] += 1
        
    class CustomProcess(QEProcess):
        def __init__(self, dir, spec, parent):
            super(QEProcess, self).__init__(dir, spec)
            self.parent = parent
        def get_results(self, analysis):
            outcome = parent.get_results(analysis)
            parent.ks.append(parent.k)
            parent.values.append(outcome)
            return outcome
