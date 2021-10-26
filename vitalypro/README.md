**Prediction of Electronic Band-structure in SiGe Superlattices**

**Overview**

The idea of this task is to train an ML model to predict electronic band-structure of a SiGe superlattice with any kind of interfacial disorder, external strain, and composition*.   
(**please see the list of assumptions that will be used in this work*)

Superlattices are periodic structures that contain layers of different materials. The interfaces significantly impact the electronic and phonon transport and make superlattices ideal candidates as thermoelectric materials. The efficiency or figure of merit of thermoelectric material is calculated as:

*zT=(S^2/ρk)T*,

where *zT* is the dimensionless figure of merit, S, ρ, k  are the Seebeck coefficient, electrical resistivity, and thermal conductivity. In order to improve the efficiency of thermoelectric material, we need to increase the Seebeck coefficient and electrical conductivity and at the same time decrease the thermal conductivity.  The former two parameters can be easily and quickly computed from the energy bands using Boltzmann equations (see https://arxiv.org/pdf/cond-mat/0602203.pdf eq. 12-16). The most time-consuming part in the computation of S is the bandstructure calculation that we will speed up with ML modeling. 
 

**Assumptions**

This task is very interesting, but it’s also complicated, and it would take a significant amount of time to complete. For this reason, I will make some assumptions that would simplify this problem:
1. I will only study superlattices that have a constant number of atoms to avoid band structure folding. As example, the electronic band structure unfolding can be done with GPAW (https://wiki.fysik.dtu.dk/gpaw/tutorialsexercises/electronic/unfold/unfold.html)
2. I will use a smaller number of atoms in the cell to speed up calculations 
3. For the transport we only care about the bands that are close to the Fermi level. Therefore, I will only predict those bands. 
4. In this work, I will only consider superlattices with the ideal interfaces. It would be interesting to investigate disordered interfaces and introduce some defects far away from the interface as well in the future work. This would require studying larger cells and introduce Voronoi tessellations to correctly describe the neighboring atoms. In this work, however, I will study the effect of external strain that comes from growing a superlattice on different substrates.
5. I will predict bandstructure along a short path (e.g. Gamma -Z) 
6. I will use the PBE exchange-correlation functional for the DFT calculations. The band gap will be underestimated but the bands shape should not be affected significantly. 

**Project details**

1. **Step 1: DFT calculations** 
....
Compute bandstructures for n SiGe SLs to build a training set for a supervised ML model. 
....

2. **Step 2: Develop a module for ML modeling**
 ...
 
3. **Step 3: Predict a SL bandstructure with the ML model**
...

