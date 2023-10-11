import functools
import urllib.request
from typing import Any, Dict, Optional, Tuple, Union

import ase.build
import ase.constraints
import numpy as np
import pymatgen.core.surface
import pymatgen.io.ase
import pymatgen.symmetry.analyzer
from exabyte_api_client.endpoints.jobs import JobEndpoints


def download_file_by_name(job_id: str, job_endpoint: JobEndpoints, target: str, pattern: str) -> None:
    """
    Downloads a file from S3 and writes it to the local disk.

    Args:
        job_id (str): The ID string for the job that the file is to be downloaded from.
        job_endpoint (JobEndpoints): Job endpoint object from the Exabyte API Client
        target (str): Target filename that the file will be written to.
        pattern (str): The name of the file as it appears on the S3. For example, to download a file named
                       "CONTCAR", the pattern would be provided as "CONTCAR"
    Notes:
        This function does not check whether the file specified in `target` exists or not, and will always
        be overwritten it if the file is successfully downloaded.
    """
    # Get a list of files for each
    job_files = job_endpoint.list_files(job_id)

    # Find the file
    for file in job_files:
        if file["name"] == pattern:
            file_metadata = file

    # Get a download URL for the file
    file_signed_url = file_metadata["signedUrl"]

    # Download the file to memory
    server_response = urllib.request.urlopen(file_signed_url)

    # Write it to disk
    with open(target, "wb") as outp:
        outp.write(server_response.read())


download_contcar = functools.partial(download_file_by_name, pattern="CONTCAR")


# Pymatgen


def is_symmetric(slab: pymatgen.core.structure.Structure) -> bool:
    """
    Checks whether a slab is in a point group with inversion symmetry, which includes the following groups:
    -1, 2/m, mmm, 4/m, 4/mmm,-3, -3m, 6/m, 6/mmm, m-3, m-3m
    This technique is borrowed from the `nonstoichiometric_symmetrized_slab` method in PyMatGen's SlabGenerator.

    Args:
        slab (pymatgen.core.structure.Structure): Slab of interest

    Returns:
        True if the slab's spacegroup has inversion symmetry, otherwise False.
    """
    # Create a Spacegroup analyzer object with the slab
    spacegroup = pymatgen.symmetry.analyzer.SpacegroupAnalyzer(slab)
    # Check for inversion symmetry
    is_symmetric = spacegroup.is_laue()
    return is_symmetric


def get_all_slabs_and_terms(
    crystal: pymatgen.core.structure.Structure, thickness: Union[int, float], is_by_layers: bool
) -> Dict[str, Dict[str, Dict[str, Any]]]:
    """
    Gets all slabs and terminations for a given crystal, forcing a specific number of layers in the resultant slab.

    Args:
        crystal (pymatgen.core.structure.Structure):  Crystal of interest
        thickness (int or float): How thick the slab is supposed to be, by either angstroms or number of layers
        is_by_layers (bool): Whether thickness is by the number of layers or by angstroms

    Returns:
        Dict
    """
    # First, get a list of all non-symetrically-equivalent indices that the crystal has
    all_indices = pymatgen.core.surface.get_symmetrically_distinct_miller_indices(crystal, max_index=3)
    slabs = {}

    for plane in all_indices:
        slab_generator = pymatgen.core.surface.SlabGenerator(
            crystal, plane, min_slab_size=thickness, min_vacuum_size=10, center_slab=True, in_unit_planes=is_by_layers
        )

        # Generate all surface terminations, and add them to the list of slabs returned
        all_terminations = slab_generator.get_slabs()
        term_dict = {}
        symmetric_terminations = filter(is_symmetric, all_terminations)
        for term, surface in enumerate(symmetric_terminations):
            ase_surface = pymatgen.io.ase.AseAtomsAdaptor.get_atoms(surface)
            term_dict[str(term)] = {"slab": ase_surface}
        slabs["".join(map(str, plane))] = term_dict
    return slabs


# ASE


def get_bulk_bottom_and_top_frac_coords(slab: ase.Atoms, layers: int = 3) -> Tuple[float, float]:
    """
    Finds the top and bottom of the bulk, in fractional coordinates

    Args:
        slab (ase.Atoms): The slab of interest
        layers (int): How many layers are in the slab

    Returns:
        A tuple containing, in order, the fractional coordintes of the bottom and top of the bulk portion of the slab
    """
    # Determine how large the slab actually is
    c_direction_coords = [atom.scaled_position[2] for atom in slab]
    slab_range = max(c_direction_coords) - min(c_direction_coords)

    # Do some math to figure out where the bottom of the bulk starts
    layer_size = slab_range / layers
    bulk_bottom = min(c_direction_coords) + layer_size

    # Figure out where the top of the bulk starts
    bulk_top = min(c_direction_coords) + 2 * layer_size

    return (bulk_bottom, bulk_top)


def freeze_center_bulk(slab: ase.Atoms) -> None:
    """
    Applies an ASE FixAtoms constraint to atoms found in the center of the slab, in-place.

    Args:
        slab (ase.Atoms): The slab of interest

    Returns:
        None, this function changes the slab in-place.
    """
    # Get the fractional coordinates for the bottom and top of the bulk
    bulk_bottom, bulk_top = get_bulk_bottom_and_top_frac_coords(slab)

    # Filter to get the atoms between the bottom/top of the bulk - that is, the bulk atoms
    frozen_atoms = filter(lambda atom: bulk_bottom <= atom.scaled_position[2] <= bulk_top, slab)

    # Get the indices of the bulk atoms, and apply the constraint to the Atoms object
    frozen_atoms_indices = [atom.index for atom in frozen_atoms]
    fix_atoms_constraint = ase.constraints.FixAtoms(indices=frozen_atoms_indices)
    slab.set_constraint(fix_atoms_constraint)


def get_vasp_total_energy(job_id: str, jobs_endpoint: JobEndpoints) -> Optional[float]:
    """
    This function takes in a VASP Job ID, reads the OUTCAR, and returns the final energy reported in the run.

    Args:
        job_id (str): The ID of the job of interest
        jobs_endpoint (JobEndpoints): A job endpoint for interacting with the Exabyte platform

    Returns:
        The electronic energy reported by the VASP job or None
    """
    # Get the URL for the OUTCAR
    files = jobs_endpoint.list_files(job_id)
    file_metadata = None
    for file in files:
        if file["name"] == "OUTCAR":
            file_metadata = file

    # Get a download URL for each CONTCAR
    cell_outcar_signed_url = file_metadata["signedUrl"]  # type: ignore

    # Download the outcar to memory
    cell_response = urllib.request.urlopen(cell_outcar_signed_url)

    # And iterate through it, finding the last electronic energy reported
    outcar = cell_response.read().decode("utf-8")
    outcar = outcar.split("\n")
    unit_cell_energy = None
    for line in outcar:
        if "sigma->0" in line:
            unit_cell_energy = float(line.strip().split()[-1])
    return unit_cell_energy


def get_surface_energy(e_slab: float, e_bulk: float, n_slab: float, n_bulk: float, a: float) -> float:
    """
    Calculates the slab energy according to the following formula:
    (E_Slab - E_bulk * (N_Slab / N_Bulk)) / (2A)
    """
    surface_energy = (e_slab - e_bulk * (n_slab / n_bulk)) / (2 * a)
    return surface_energy


def get_slab_area(a_vector: np.ndarray, b_vector: np.ndarray) -> float:
    """
    Gets the area of a slab defined by the two unit vectors.
    The magnitude of the cross product of the vectors is the area of the parallelogram they enclose.

    Args:
        a_vector
    """
    crossprod = np.cross(a_vector, b_vector)
    magnitude = np.linalg.norm(crossprod)
    return magnitude
