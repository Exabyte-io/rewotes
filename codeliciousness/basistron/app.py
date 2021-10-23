
# -*- coding: utf-8 -*-
import sys
from typing import List

import pandas as pd

from basistron import basis, cccbdb, cli, exabyte, utils

log = utils.get_logger("basistron.app")


def filter_dfs_by_name(
    dfs: List[pd.DataFrame],
    regime: str,
    match: str = "standard",
):
    """Localize heuristic table selection here."""
    log.info(f"choosing dataframe from {[df.name for df in dfs]}")
    this = None
    if len(dfs) == 1:
        this = dfs[0]
    elif regime == "calculated":
        # calculated results usually show up with
        # three tables of result groups
        # ["empirical", "standard", "effective"]
        for df in dfs:
                if match in df.name:
                    if this is not None:
                        raise Exception("duplicate match logic")
                    this = df
                    break
    if this is None:
        raise Exception("could not find reference data")
    return this


def main(args):
    """Run the BasisTron 5000!"""

    # init
    driver = cli.process_args(cli.get_parser().parse_args(args))
    formula = driver.simple_formula()
    log.info(f"starting basis set selector on {formula} for {driver.property}")

    # refdata
    db = cccbdb.Cccbdb()
    dfs = db.get_dataframes(formula, driver.property.value)
    if not dfs:
        log.error("found no tables from CCCBDB")
        sys.exit()
    try:
        df = filter_dfs_by_name(dfs, driver.regime.value)
    except Exception as e:
        log.error("failed to select table from reference data")
        sys.exit()

    # select reference datum
    if driver.value is None:
        print(df)
        allowed_basis_sets = []
    else:
        # sanity check provided reference data
        lower, upper = driver.acceptable_range()
        acceptable = df[(df >= lower) & (df <= upper)]
        if not acceptable.sum().sum():
            msg = f"{driver.value}±{driver.tolerance:.2f}%"
            log.error(f"no reference data found within {msg}")
            log.warning("subsequent analysis may fail")
        allowed_basis_sets = []

    # basis set analysis
    # Of allowed_basis_sets, pick the most compact
    bases = basis.Basis.load_basis_sets()
    ranked = basis.Basis.rank_basis_sets(bases)
    total_allowed = basis.Basis.get_allowed_basis_sets(
        ranked, list(set([r[0] for r in driver.xyz_data]))
    )
    log.info(f"found {len(allowed_basis_sets)} allowed basis sets")
    log.info(f"ranking allowed out of {len(total_allowed)} basis sets")

    # update workflow with basis set, level of theory?

    # ship it
    ebc = exabyte.Client()
    config = ebc.get_material_config("", driver.xyz_data_to_dict())
    material = ebc.get_endpoint("material").create(config)
    workflow = ebc.get_workflow()
    job_cfg = ebc.get_job_config(
        workflow["owner"]["_id"],
        material["_id"],
        workflow["_id"],
        "basistron.app",
    )
    print(job_cfg)
    # job = c.submit_job(job_cfg)


if __name__ == "__main__":
    main(sys.argv[1:])