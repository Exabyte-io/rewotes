# -*- coding: utf-8 -*-
import json
import logging
import sys
from typing import Any, Dict, List, Optional

import pandas as pd

from basistron import basis, cccbdb, cli, exabyte, model, utils

log = utils.get_logger("basistron.app")


class MissingReferenceData(Exception):
    pass


def filter_dfs_by_name(
    dfs: List[pd.DataFrame],
    regime: str,
    match: str = "standard",
) -> pd.DataFrame:
    """Localize heuristic table selection here."""
    log.info(f"choosing dataframe from {[df.name[:20] for df in dfs]}")
    this = None
    if len(dfs) == 1:
        this = dfs[0]
    # calculated results usually show up with
    # three tables of result groups
    # ["empirical", "standard", "effective"]
    elif regime == "calculated":
        for df in dfs:
            if match in df.name:
                if this is not None:
                    raise Exception("duplicate match logic")
                this = df
                break
    # table structure very different for exptl data from CCCBDB
    # so would need its own logic presumably
    # elif regime == "experimental"
    if this is None:
        raise MissingReferenceData("table filtering logic failed")
    return this


def set_reference_value(
    df: pd.DataFrame,
    driver: model.Execution,
) -> None:
    """Assume input data is ordered approximately as follows:
    Increasing index (ordered) values means increasing level of theory.
    Increasing column (ordered) values means increasing basis quality.
    """
    theory = driver.reference_theory
    basis = driver.reference_basis
    log.info(f"looking for reference datum @ ({theory},{basis})")
    if theory not in df.index:
        log.warning(f"reference theory {theory} not in {df.index.values}")
        # should probably pick "best" reference theory with known available data
        theory = df.index.values[-1]
        log.info(f"selecting best reference theory available: {theory}")
    # case insensitive comparisons but this could be better
    if basis not in df.columns:
        log.warning(f"reference basis {basis} not found in {df.columns.values}")
        # CCCBDB order is not strictly increasing so this is a hack
        basis = df.columns[len(df.columns) // 2 - 2]
        log.info(f"selecting medium size reference basis set: {basis}")
    value = df.loc[theory, basis]
    if pd.isnull(value):
        log.error(f"reference datum for ({theory},{basis}) not found")
        raise MissingReferenceData("reference datum selection logic failed")
    log.info(f"selected reference value = {value} @ ({theory},{basis})")
    driver.reference_theory = theory
    driver.reference_basis = basis
    driver.value = value


def select_basis_set(
    df: pd.DataFrame,
    driver: model.Execution,
) -> str:
    """Filter allowed basis sets for to determine the best available.
    If a basis set cannot be found for the target level of thery,
    a target level of theory will be chosen if possible. In the case
    of multiple available basis sets, attempt to choose the most compact
    one from the basis set database. If unavailable, assume the basis
    sets are already ordered in increasing size.
    """

    def get_basis_sets(
        df: pd.DataFrame, theory: str, lower: float, upper: float
    ) -> Optional[str]:
        try:
            target = df.loc[theory]
        except KeyError:
            return
        acceptable = target[(target >= lower) & (target <= upper)]
        if acceptable.any():
            log.info(f"selecting {len(acceptable)} basis sets for {theory}")
            log.info(f"bounds on selected ({acceptable.min()},{acceptable.max()})")
            return acceptable.index.values.tolist()

    allowed = f"within {driver.value}±{driver.tolerance:.2f}%"
    lower, upper = driver.acceptable_range()
    # allow priority to user provided target theory
    theories = [driver.target_theory]
    theories = theories + df.index.difference(theories).tolist()
    for theory in theories:
        log.info(f"attempting to select basis set for target theory {theory}")
        basis_sets = get_basis_sets(df, theory, lower, upper)
        if basis_sets is not None:
            driver.target_theory = theory
            break

    if basis_sets is None:
        msg = f"could not find any basis set for any theory to yield results {allowed}"
        log.error(msg)
        raise MissingReferenceData(f"no target data found {allowed}")

    try:
        # basis set analysis doomed to fail for now
        bases = basis.Basis.load_basis_sets()
        ranked = basis.Basis.rank_basis_sets(bases)
        # needs more cleanup to match off against CCCBDB basis set specs
        total_allowed = [
            ".".join(b.lower().split(".")[:-1])
            for b in basis.Basis.get_allowed_basis_sets(
                ranked, list(set([r[0] for r in driver.xyz_data]))
            )
        ]
        ordered = [basis for basis in total_allowed if basis in basis_sets]
    except Exception:
        ordered = []

    if not ordered:
        return basis_sets[0]
    log.warning("did not find matching basis sets in database, skipping")
    return ordered[0]


def main(args):
    """Run the BasisTron 5000! Business logic is broken into four
    main parts.

        1. Filter the CCCBDB data tables that are provided.
           It is use-case specific enough to belong here.
        2. Determine the reference value (if not provided),
           keeping track of the level of theory and basis
        3. Choose a basis set at a target level of theory
           providing the same accuracy within the reference
           tolerance.
        4. Create the material, update the workflow, submit
           the job to the exabyte cluster.
    """

    driver = cli.process_args(cli.get_parser().parse_args(args))
    formula = driver.simple_formula()
    log.info(f"starting basis set selector on {formula} for {driver.property}")

    # step 1
    db = cccbdb.Cccbdb()
    dfs = db.get_dataframes(formula, driver.property.value)
    if not dfs:
        log.error("found no tables from CCCBDB")
        return
    try:
        df = filter_dfs_by_name(dfs, driver.regime.value)
    except MissingReferenceData as e:
        log.error(f"failed to select table from reference data: {repr(e)}")
        return

    # step 2
    if driver.value is None:
        set_reference_value(df, driver)
    else:
        # reference_theory and reference_basis are unused
        # sanity check provided reference value
        lower, upper = driver.acceptable_range()
        acceptable = df[(df >= lower) & (df <= upper)]
        if not acceptable.any().any():
            msg = f"{driver.value}±{driver.tolerance:.2f}%"
            log.error(f"no reference data found within {msg}")
            log.warning("subsequent analysis may fail")

    # step 3
    orig = driver.target_theory
    basis_set = select_basis_set(df, driver)
    curr = driver.target_theory
    if orig != curr:
        msg = f"updated target theory from {orig} to {curr} to meet tolerance"
        log.warning(msg)

    # step 4
    def debug(name: str, blob: Dict[str, Any]):
        for ln in json.dumps(config, indent=4).splitlines():
            log.debug(f"{name}: {ln}")

    log.setLevel(logging.DEBUG)

    ebc = exabyte.Client()
    config = ebc.get_material_config("", driver.xyz_data_to_dict())
    debug("config", config)
    material = ebc.get_endpoint("material").create(config)
    debug("material", material)
    workflow = ebc.get_workflow()
    debug("workflow", workflow)
    job_cfg = ebc.get_job_config(
        workflow["owner"]["_id"],
        material["_id"],
        workflow["_id"],
        "basistron.app",
    )
    debug("job config", job_cfg)
    job = ebc.submit_job(job_cfg)
    debug("job", job)
    # print(json.dumps(job, indent=4))
    # log.info(f"successfully submitted job with ID={job['_id']}")


if __name__ == "__main__":
    main(sys.argv[1:])
