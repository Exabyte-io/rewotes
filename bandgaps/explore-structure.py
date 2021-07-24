import os
import numpy as np

import matplotlib.pyplot as plt

from pymatgen.core.structure import Structure
from pymatgen.transformations.advanced_transformations import EnumerateStructureTransformation

from megnet.utils.models import load_model


PATH = os.path.dirname(os.path.abspath(__file__))

# odst = SQSTransformation()
order = EnumerateStructureTransformation(max_cell_size=10, sort_criteria="nsites")

si_sc = Structure.from_file(PATH + f"/data/POSCAR.mp-149_Si")

alloys = []
for x in np.linspace(0.05, 0.95, 18):
    tmp = si_sc.copy()
    tmp.replace_species({"Si": {"Si": 1 - x, "Ge": x}})
    tmp = order.apply_transformation(tmp)

    alloys.append(tmp)

pt_model = load_model("Bandgap_MP_2018")

y_pred = pt_model.predict_structures(alloys).ravel()


plt.figure(figsize=(8, 8))
plt.xlim
plt.plot(np.linspace(0.05, 0.95, 18), y_pred, label="MEGNet Prediction")

# https://journals.aps.org/prb/abstract/10.1103/PhysRevB.40.5683
fitted = [1.155 - 0.43*x + 0.206*x**2 if x < 0.85 else 2.010 - 1.270*x for x in np.linspace(0, 1)]

plt.plot(np.linspace(0, 1), fitted, label="Experimental Fit")
ax = plt.gca()
ax.set_xlim(0, 1)
ax.set_ylim(0, 1.5)

# plt.scatter((0, 0.125, 0.25, 0.5, 0.75, 1), (0.8527, 0.8585, 0.0000, 0.4267, 0.0000, 0.0000), color="red", marker="x", clip_on=False, zorder=3, label="Training Data")
# plt.scatter((0, 0.125, 0.5), (0.8527, 0.8585, 0.4267), color="red", marker="x", clip_on=False, zorder=3, label="Training Data")

plt.legend(frameon=False)

xlabel = r"$x$" + " in " + r"$Si_{x}Ge_{1-x}$"
ax.set_xlabel(xlabel)

ylabel = r"$E_g$" + " / eV"
ax.set_ylabel(ylabel)

plt.show()