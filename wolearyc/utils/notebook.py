import os
import re
from urllib.parse import unquote, urlparse

import requests


def get_notebook_info() -> dict:
    """
    Get the information about a currently running notebook in Google Colab.
    Args:
        None
    Return:
        a dict with notebook info.
    """
    # ip = socket.gethostbyname(socket.gethostname())  # 172.28.0.12
    ip = os.getenv("COLAB_JUPYTER_IP")
    response = requests.get(f"http://{ip}:9000/api/sessions").json()[0]

    notebook_name = response["name"]
    path = urlparse(unquote(response["path"]).split("=")[1]).path
    parsed = re.findall("(.*)/blob/(.*)/examples/(.*)", path)[0]
    github_org_repo = parsed[0][1:]  # remove leading /
    branch_name = parsed[1]
    notebook_path = f"examples/{parsed[2]}"

    return dict(
        notebook_name=notebook_name,
        notebook_path=notebook_path,
        branch_name=branch_name,
        github_org_repo=github_org_repo,
    )


def print_notebook_path() -> None:
    """
    A proxy function used for a single string return when
    the corresponding entry-point script in 'setup.py' is called.

    Args:
        None
    Return:
        None
    """
    # 'return get_notebook_info()["notebook_path"]' returns a non-zero exit code; using 'print' instead.
    print(get_notebook_info()["notebook_path"])
