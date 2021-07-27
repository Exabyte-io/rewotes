# %%
import os
import numpy as np

import matplotlib.pyplot as plt

from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.transformations.advanced_transformations import (
    EnumerateStructureTransformation,
)

from tqdm import tqdm
from megnet.utils.models import load_model


PATH = os.path.dirname(os.path.abspath(__file__))

# %%

order = EnumerateStructureTransformation(
    max_cell_size=10,
    sort_criteria="nsites",
)

si_sc = Structure.from_file(PATH + f"/data/POSCAR.mp-149_Si")
analyzer = SpacegroupAnalyzer(si_sc)
si_sc = analyzer.find_primitive()

# %%
def swap_si_ge(struct):
    """swap Si and Ge in structure

    Args:
        struct (Structure): a pymatgen structure

    Returns:
        structure with Si and Ge swapped
    """
    swap = {"Si": {"Ge": 1}, "Ge": {"Si": 1}}

    tmp = struct.copy()
    tmp.replace_species(swap)

    return tmp


alloys = [si_sc, swap_si_ge(si_sc)]
for x in tqdm([0.1, 0.2, 0.3, 0.4, 0.5]):  # takes ~ 15 minutes to run
    tmp = si_sc.copy()
    tmp.replace_species({"Si": {"Si": 1 - x, "Ge": x}})
    # NOTE we get an error in the subprocess execution sadly the error message is uninformative
    # further investigation needed.
    tmp = order.apply_transformation(tmp, return_ranked_list=10)  # limit to 10 for evaluation purposes

    alloys.extend(s["structure"] for s in tmp)

    if x != 0.5:
        ge_tmp = [swap_si_ge(s["structure"]) for s in tmp]
        alloys.extend(ge_tmp)

# %%

pt_model = load_model("Bandgap_MP_2018")

# %%

x_si = [1 - x.composition.get_atomic_fraction("Si") for x in alloys]

# [s.scale_lattice(v) for v, s in zip(vol_alloy, alloys)]

y_pred = pt_model.predict_structures(alloys).ravel()

# TODO load energy based pre_trained model, take boltzman weighted average of band gap, plot the uncertainty of
# the boltzman ensemble.

plt.figure(figsize=(8, 8))
plt.xlim
plt.scatter(x_si, y_pred, label="MEGNet Prediction", alpha=0.3, clip_on=False, zorder=3)

# https://journals.aps.org/prb/abstract/10.1103/PhysRevB.40.5683
fitted = [1.155 - 0.43*x + 0.206*x**2 if x < 0.85 else 2.010 - 1.270*x for x in np.linspace(0, 1)]

plt.plot(np.linspace(0, 1), fitted, color="tab:orange", label="Experimental Fit")
ax = plt.gca()
ax.set_xlim(0, 1)
ax.set_ylim(0, 2.5)

plt.legend(frameon=False)

xlabel = r"$x$" + " in " + r"$Si_{1-x}Ge_{x}$"
ax.set_xlabel(xlabel)

ylabel = r"$E_g$" + " / eV"
ax.set_ylabel(ylabel)
plt.savefig(PATH+f"/SiGe-disorded-structures.pdf")
plt.show()
# %%

# TODO scale the volume of the crystal
# Si = 3.867 \AA, 40.888 - mp-149
# SiGe = 3.955 \AA, 43.741 - mp-1219182
# Ge = 3.955 \AA, 47.847 - mp-32

def get_est_vol(x):
    """[summary]

    Args:
        x (float): fraction of Ge in structure

    Returns:
        estimated volume of unit cell
    """
    return (40.888 * (1 - x) + 43.741 * x) if x < 0.5 else (43.741 * (1 - x) + 47.847 * x)


vol_alloy = [get_est_vol(x) for x in x_si]

[s.scale_lattice(v) for v, s in zip(vol_alloy, alloys)]

y_pred = pt_model.predict_structures(alloys).ravel()

plt.figure(figsize=(8, 8))
plt.xlim
plt.scatter(x_si, y_pred, label="MEGNet Prediction", alpha=0.3, clip_on=False, zorder=3)

# https://journals.aps.org/prb/abstract/10.1103/PhysRevB.40.5683
fitted = [1.155 - 0.43*x + 0.206*x**2 if x < 0.85 else 2.010 - 1.270*x for x in np.linspace(0, 1)]

plt.plot(np.linspace(0, 1), fitted, color="tab:orange", label="Experimental Fit")
ax = plt.gca()
ax.set_xlim(0, 1)
ax.set_ylim(0, 2.5)

plt.legend(frameon=False)

xlabel = r"$x$" + " in " + r"$Si_{1-x}Ge_{x}$"
ax.set_xlabel(xlabel)

ylabel = r"$E_g$" + " / eV"
ax.set_ylabel(ylabel)
plt.savefig(PATH+f"/SiGe-disorded-structures-scaled.pdf")
plt.show()