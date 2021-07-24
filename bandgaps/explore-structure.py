
import os
import numpy as np

import matplotlib.pyplot as plt

from pymatgen.core.structure import Structure
from pymatgen.transformations.advanced_transformations import EnumerateStructureTransformation

from megnet.utils.models import load_model


PATH = os.path.dirname(os.path.abspath(__file__))


order = EnumerateStructureTransformation(max_cell_size=10, sort_criteria="nsites")

si_sc = Structure.from_file(PATH + f"/data/POSCAR.mp-149_Si")

alloys = []
for x in np.linspace(0.05, 0.95, 18):
    tmp = si_sc.copy()
    tmp.replace_species({"Si": {"Si": 1 - x, "Ge": x}})
    # NOTE we get an error in the subprocess execution sadly the error message is uninformative
    # further investigation needed.
    tmp = order.apply_transformation(tmp)

    alloys.append(tmp)

pt_model = load_model("Bandgap_MP_2018")

y_pred = pt_model.predict_structures(alloys).ravel()

# TODO load energy based pre_trained model, take boltzman weighted average of band gap, plot the uncertainty of
# the boltzman ensemble.

