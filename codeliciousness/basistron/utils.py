# -*- coding: utf-8 -*-
import os
import logging

logging.basicConfig()


def get_logger(name, level=logging.INFO):
    log = logging.getLogger(name)
    log.setLevel(level)
    return log


def default_cache_dir():
    path = os.path.join(
        os.path.expanduser("~"),
        ".basistron",
    )
    os.makedirs(path, exist_ok=True)
    return path


class Log:

    @property
    def log(self):
        return get_logger(
            ".".join([
                self.__module__,
                self.__class__.__name__,
            ])
        )

class _env:
    """Namespace collecting all environment variables
    used within the application."""

    @property
    def exabyte_host(self):
        return os.getenv("EXABYTE_HOST", "platform.exabyte.io")

    @property
    def exabyte_port(self):
        return os.getenv("EXABYTE_PORT", 443)

    @property
    def exabyte_username(self):
        return os.getenv("EXABYTE_USERNAME")

    @property
    def exabyte_password(self):
        return os.getenv("EXABYTE_PASSWORD")

    @property
    def exabyte_client_id(self):
        """Used as X-Account-Id header"""
        return os.getenv("EXABYTE_CLIENT_ID")

    @property
    def exabyte_client_secret(self):
        """Used as X-Account-Id header"""
        return os.getenv("EXABYTE_CLIENT_SECRET")

env = _env()
