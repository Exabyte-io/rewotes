from . import seekpath
from markupsafe import escape
import copy
import io
import json
import time
import traceback

import numpy as np

from ase.data import chemical_symbols
from tools_barebone import logme, get_tools_barebone_version
from tools_barebone.structure_importers import get_structure_tuple, UnknownFormatError




MAX_NUMBER_OF_ATOMS = 1000


def get_json_for_visualizer(
    cell, relcoords, atomic_numbers
):  # pylint: disable=too-many-locals
    from seekpath import hpkot, brillouinzone
    hpkot = seekpath.hpkot

    system = (np.array(cell), np.array(relcoords), np.array(atomic_numbers))
    res = hpkot.get_path(system, with_time_reversal=False)

    real_lattice = res["primitive_lattice"]
    # rec_lattice = np.linalg.inv(real_lattice).T # Missing 2pi!
    rec_lattice = np.array(hpkot.tools.get_reciprocal_cell_rows(real_lattice))
    b1, b2, b3 = rec_lattice
    faces_data = seekpath.brillouinzone.get_BZ(b1=b1, b2=b2, b3=b3)

    response = {}
    response["faces_data"] = faces_data
    response["b1"] = b1.tolist()
    response["b2"] = b2.tolist()
    response["b3"] = b3.tolist()
    ## Convert to absolute
    response["kpoints"] = {
        k: (v[0] * b1 + v[1] * b2 + v[2] * b3).tolist()
        for k, v in res["point_coords"].items()
    }
    response["kpoints_rel"] = {
        k: [v[0], v[1], v[2]] for k, v in res["point_coords"].items()
    }
    response["path"] = res["path"]

    # It should use the same logic, so give the same cell as above
    res_explicit = seekpath.get_explicit_k_path(system, with_time_reversal=False)
    for k in res_explicit:
        if k == "segments" or k.startswith("explicit_"):
            if isinstance(res_explicit[k], np.ndarray):
                response[k] = res_explicit[k].tolist()
            else:
                response[k] = res_explicit[k]

    if (
        np.sum(
            np.abs(
                np.array(res_explicit["reciprocal_primitive_lattice"])
                - np.array(res["reciprocal_primitive_lattice"])
            )
        )
        > 1.0e-7
    ):
        raise AssertionError("Got different reciprocal cells...")

    # Response for JS, and path_results
    return response, res


def process_structure_core(  # pylint: disable=too-many-locals,too-many-statements,too-many-arguments
    filecontent,
    fileformat,
    call_source="",
    logger=None,
    flask_request=None,
):
    
    start_time = time.time()
    fileobject = io.StringIO(str(filecontent))
    # form_data = dict(flask_request.form)
    form_data = None
    try:
        structure_tuple = get_structure_tuple(
            fileobject, fileformat, extra_data=None
        )
    except UnknownFormatError:
        logme(
            logger,
            filecontent,
            fileformat,
            flask_request,
            call_source,
            reason="unknownformat",
            extra={"form_data": form_data,},
        )
        raise FlaskRedirectException("Unknown format '{}'".format(fileformat))
    except Exception:
        # There was an exception...
        logme(
            logger,
            filecontent,
            fileformat,
            flask_request,
            call_source,
            reason="exception",
            extra={"traceback": traceback.format_exc(), "form_data": form_data,},
        )
        raise FlaskRedirectException(
            "I tried my best, but I wasn't able to load your "
            "file in format '{}'...".format(fileformat)
        )

    if len(structure_tuple[1]) > MAX_NUMBER_OF_ATOMS:
        ## Structure too big
        logme(
            logger,
            filecontent,
            fileformat,
            flask_request,
            call_source,
            reason="toolarge",
            extra={"number_of_atoms": len(structure_tuple[1]), "form_data": form_data,},
        )
        raise FlaskRedirectException(
            "Sorry, this online visualizer is limited to {} atoms "
            "in the input cell, while your structure has {} atoms."
            "".format(MAX_NUMBER_OF_ATOMS, len(structure_tuple[1]))
        )

    # Log the content in case of valid structure
    logme(
        logger,
        filecontent,
        fileformat,
        flask_request,
        call_source,
        reason="OK",
        extra={"number_of_atoms": len(structure_tuple[1]), "form_data": form_data,},
    )

    try:
        in_json_data = {
            "cell": structure_tuple[0],
            "scaled_coords": structure_tuple[1],
            "atomic_numbers": structure_tuple[2],
        }

        out_json_data, path_results = get_json_for_visualizer(
            in_json_data["cell"],
            in_json_data["scaled_coords"],
            in_json_data["atomic_numbers"],
        )

        raw_code_dict = copy.copy(out_json_data)
        for k in list(raw_code_dict.keys()):
            if k.startswith("explicit_"):
                raw_code_dict.pop(k)
            if k == "segments":
                raw_code_dict.pop(k)
        raw_code_dict.pop("faces_data")
        raw_code_dict["primitive_lattice"] = path_results["primitive_lattice"].tolist()
        raw_code_dict["primitive_positions"] = path_results[
            "primitive_positions"
        ].tolist()
        inputstructure_positions_cartesian = np.dot(
            np.array(in_json_data["scaled_coords"]), np.array(in_json_data["cell"]),
        ).tolist()
        primitive_positions_cartesian = np.dot(
            np.array(path_results["primitive_positions"]),
            np.array(path_results["primitive_lattice"]),
        ).tolist()
        primitive_positions_cartesian_refolded = np.dot(
            np.array(path_results["primitive_positions"]) % 1.0,
            np.array(path_results["primitive_lattice"]),
        ).tolist()
        raw_code_dict["primitive_positions_cartesian"] = primitive_positions_cartesian

        # raw_code['primitive_types'] = path_results['primitive_types']
        primitive_symbols = [
            chemical_symbols[num] for num in path_results["primitive_types"]
        ]
        raw_code_dict["primitive_symbols"] = primitive_symbols

        raw_code = json.dumps(raw_code_dict, indent=2)
        ## I manually escape it to then add <br> and pass it to a filter with
        ## |safe. I have to 'unicode' it otherwise it keeps escaping also the
        ## next replaces
        raw_code = (
            str(escape(raw_code)).replace("\n", "<br>").replace(" ", "&nbsp;")
        )

        kpoints = [
            [
                k,
                out_json_data["kpoints"][k][0],
                out_json_data["kpoints"][k][1],
                out_json_data["kpoints"][k][2],
            ]
            for k in sorted(out_json_data["kpoints"])
        ]
        kpoints_rel = [
            [
                k,
                out_json_data["kpoints_rel"][k][0],
                out_json_data["kpoints_rel"][k][1],
                out_json_data["kpoints_rel"][k][2],
            ]
            for k in sorted(out_json_data["kpoints_rel"])
        ]

        inputstructure_cell_vectors = [
            [idx, coords[0], coords[1], coords[2]]
            for idx, coords in enumerate(in_json_data["cell"], start=1)
        ]
        inputstructure_symbols = [
            chemical_symbols[num] for num in in_json_data["atomic_numbers"]
        ]
        inputstructure_atoms_scaled = [
            [label, coords[0], coords[1], coords[2]]
            for label, coords in zip(
                inputstructure_symbols, in_json_data["scaled_coords"]
            )
        ]
        inputstructure_atoms_cartesian = [
            [label, coords[0], coords[1], coords[2]]
            for label, coords in zip(
                inputstructure_symbols, inputstructure_positions_cartesian
            )
        ]

        direct_vectors = [
            [idx, coords[0], coords[1], coords[2]]
            for idx, coords in enumerate(path_results["primitive_lattice"], start=1)
        ]

        reciprocal_primitive_vectors = [
            [idx, coords[0], coords[1], coords[2]]
            for idx, coords in enumerate(
                path_results["reciprocal_primitive_lattice"], start=1
            )
        ]

        atoms_scaled = [
            [label, str(coords[0]), str(coords[1]), str(coords[2])]
            for label, coords in zip(
                primitive_symbols, path_results["primitive_positions"]
            )
        ]

        atoms_cartesian = [
            [label, coords[0], coords[1], coords[2]]
            for label, coords in zip(primitive_symbols, primitive_positions_cartesian)
        ]

        # Create extetically-nice looking path, with dashes and pipes
        suggested_path = []
        if path_results["path"]:
            suggested_path.append(path_results["path"][0][0])
            suggested_path.append("-")
            suggested_path.append(path_results["path"][0][1])
            last = path_results["path"][0][1]
        for p1, p2 in path_results["path"][1:]:
            if p1 != last:
                suggested_path.append("|")
                suggested_path.append(p1)
            suggested_path.append("-")
            suggested_path.append(p2)
            last = p2

        primitive_lattice = path_results["primitive_lattice"]
        xsfstructure = []
        xsfstructure.append("CRYSTAL")
        xsfstructure.append("PRIMVEC")
        for vector in primitive_lattice:
            xsfstructure.append("{} {} {}".format(vector[0], vector[1], vector[2]))
        xsfstructure.append("PRIMCOORD")
        xsfstructure.append("{} 1".format(len(primitive_positions_cartesian_refolded)))
        for atom_num, pos in zip(
            path_results["primitive_types"], primitive_positions_cartesian_refolded
        ):
            xsfstructure.append("{} {} {} {}".format(atom_num, pos[0], pos[1], pos[2]))
        xsfstructure = "\n".join(xsfstructure)

        compute_time = time.time() - start_time
    except Exception:
        logme(
            logger,
            filecontent,
            fileformat,
            flask_request,
            call_source,
            reason="codeexception",
            extra={"traceback": traceback.extract_stack(), "form_data": form_data,},
        )
        raise

    return dict(
        jsondata=json.dumps(out_json_data),
        volume_ratio_prim=int(round(path_results["volume_original_wrt_prim"])),
        raw_code=raw_code,
        kpoints=kpoints,
        kpoints_rel=kpoints_rel,
        bravais_lattice=path_results["bravais_lattice"],
        bravais_lattice_extended=path_results["bravais_lattice_extended"],
        spacegroup_number=path_results["spacegroup_number"],
        spacegroup_international=path_results["spacegroup_international"],
        direct_vectors=direct_vectors,
        inputstructure_cell_vectors=inputstructure_cell_vectors,
        inputstructure_atoms_scaled=inputstructure_atoms_scaled,
        inputstructure_atoms_cartesian=inputstructure_atoms_cartesian,
        atoms_scaled=atoms_scaled,
        with_without_time_reversal=(
            "with" if path_results["has_inversion_symmetry"] else "without"
        ),
        atoms_cartesian=atoms_cartesian,
        reciprocal_primitive_vectors=reciprocal_primitive_vectors,
        primitive_vectors=primitive_lattice,
        suggested_path=suggested_path,
      )


def primitive(fileContent,fileformat):
    with open(fileContent) as file:
        fileContent = file.read()
    process_structure_core(filecontent=fileContent,fileformat=fileformat)
    structure = process_structure_core(filecontent=fileContent,fileformat=fileformat)
    lattice = (structure["primitive_vectors"].astype('str')).tolist()
    atoms = structure["atoms_scaled"]
    kpoints = structure['kpoints_rel']
    return (lattice,atoms,kpoints)



