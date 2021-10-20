
# -*- coding: utf-8 -*-
import os
import sys

from basistron import cli
from basistron import utils
from basistron import client

log = utils.get_logger("basistron.app")


def main(args):

    parser = cli.get_parser()
    driver = cli.process_args(parser.parse_args(args))
    log.info(f"starting basis set selector for {driver.target_property}")

    c = client.Client()
    config = c.get_material_config("", driver.xyz_data_to_dict())
    material = c.get_endpoint("material").create(config)
    workflow = c.get_workflow()
    job_cfg = c.get_job_config(
        workflow["owner"]["_id"],
        material["_id"],
        workflow["_id"],
        "basistron.app",
    )
    job = c.submit_job(job_cfg)



if __name__ == "__main__":
    main(sys.argv[1:])


