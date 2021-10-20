# -*- coding: utf-8 -*-
import os
from types import ModuleType
from typing import Dict

from exabyte_api_client import endpoints

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

	_endpoints = _collect_endpoints()

	def get_endpoint(self, name: str):
		endpoint = self._endpoints.get(name)
		if endpoint is None:
			self.log.warning(f"name {name} not found in {self._endpoints.keys()}")
			return
		args = (
			utils.env.exabyte_host,
			utils.env.exabyte_port,
			utils.env.exabyte_username,
			utils.env.exabyte_password,
		) if name == "login" else (
			utils.env.exabyte_host,
			utils.env.exabyte_port,
			utils.env.exabyte_client_id,
			utils.env.exabyte_client_secret,
		)
		return endpoint(*args)

	@property
	def owner_query(self):
		return {"owner._id": utils.env.exabyte_client_id}

	@property
	def default_query(self):
		return {"isDefault": True, **self.owner_query}

	def __init__(self):
		tokens = self.get_endpoint("login").login()
		os.environ["EXABYTE_CLIENT_ID"] = tokens["X-Account-Id"]
		os.environ["EXABYTE_CLIENT_SECRET"] = tokens["X-Auth-Token"]


if __name__ == "__main__":
	c = Client()
	project = c.get_endpoint("project")
	[project_data] = project.list(c.default_query)
	project_id = project_data["_id"]
	owner_id = project_data["owner"]["_id"]
	print("default project data", project_id, owner_id)
	workflows = c.get_endpoint("workflow")
	print(len(workflows.list(c.default_query)))
