# -*- coding: utf-8 -*-
import os
from types import ModuleType
from typing import Any, Dict, List, Union

from exabyte_api_client import endpoints
from requests import HTTPError

from basistron import utils

log = utils.get_logger(__name__)


def _collect_endpoints() -> Dict[str, ModuleType]:
    """Collect all the exabyte clients and
    wrap them in a single API Client."""

    def _iter_module(module: ModuleType):
        for k, v in vars(module).items():
            if k.startswith("__"):
                continue
            if k.isupper():
                continue
            if k in ["json"]:
                continue
            yield k, v

    apis = {}
    # TODO : missing charges endpoint
    for _, val in _iter_module(endpoints):
        if isinstance(val, ModuleType):
            for sub, cls in _iter_module(val):
                if sub.startswith("Base"):
                    continue
                if sub.endswith("Endpoint") or sub.endswith("Endpoints"):
                    key = sub.replace("Endpoints", "").replace("Endpoint", "").lower()
                    apis[key] = cls
    return apis


class Client(utils.Log):
    """Wrapper around endpoints API for simplicity."""

    __endpoints = _collect_endpoints()

    def get_endpoint(self, name: str) -> endpoints.BaseEndpoint:
        """Return an instance of an exabyte Endpoint and store
        it for subsequent calls to the same endpoint."""
        endpoint = self._endpoints.get(name)
        if endpoint is not None:
            return endpoint
        endpoint = self.__endpoints.get(name)
        if endpoint is None:
            self.log.warning(f"name {name} not found in {self._endpoints.keys()}")
            return
        args = (
            (
                utils.env.exabyte_host,
                utils.env.exabyte_port,
                utils.env.exabyte_username,
                utils.env.exabyte_password,
            )
            if name == "login"
            else (
                utils.env.exabyte_host,
                utils.env.exabyte_port,
                utils.env.exabyte_client_id,
                utils.env.exabyte_client_secret,
            )
        )
        self._endpoints[name] = endpoint(*args)
        return self._endpoints[name]

    def get_project(self, query: Dict[str, Any] = None) -> Dict[str, Any]:
        project = self.get_endpoint("project")
        [project_data] = project.list(query or self.default_query)
        return project_data

    def get_workflow(self, query: Dict[str, Any] = None) -> Dict[str, Any]:
        workflow = self.get_endpoint("workflow")
        [workflow_data] = workflow.list(query or self.default_query)
        return workflow_data

    def get_input_template(self) -> Dict[str, Any]:
        workflow = self.get_workflow()
        [unit] = workflow["subworkflows"][0]["units"]
        return unit["input"][0]

    def get_job_config(
        self,
        owner_id: str,
        material_id: str,
        workflow_id: str,
        job_name: str,
    ) -> Dict[str, Union[str, Dict[str, str]]]:
        return {
            "owner": {
                "_id": owner_id,
            },
            "_material": {
                "_id": material_id,
            },
            "workflow": {
                "_id": workflow_id,
            },
            "name": job_name,
        }

    def get_material_config(
        self,
        name: str,
        basis: Dict[str, List[Dict[str, Any]]],
    ) -> Dict[str, Any]:
        return {
            "name": name,
            "basis": {
                "units": "cartesian",
                "name": "basis",
                **basis,
            },
            "tags": ["basistron"],
        }

    def submit_job(self, config: Dict[str, str]):
        jobs = self.get_endpoint("job")
        job = jobs.create(config)
        jobs.submit(job["_id"])
        return job

    @property
    def owner_query(self) -> Dict[str, str]:
        return {"owner._id": utils.env.exabyte_client_id}

    @property
    def default_query(self) -> Dict[str, Any]:
        return {"isDefault": True, **self.owner_query}

    def __init__(self):
        self._endpoints = {}
        try:
            tokens = self.get_endpoint("login").login()
            os.environ["EXABYTE_CLIENT_ID"] = tokens["X-Account-Id"]
            os.environ["EXABYTE_CLIENT_SECRET"] = tokens["X-Auth-Token"]
        except HTTPError:
            self.log.error(
                "authentication failure, client is useless. are "
                "EXABYTE_USERNAME+EXABYTE_PASSWORD in the environment?"
            )
