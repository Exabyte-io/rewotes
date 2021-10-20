
from basistron import client, utils


def test_client(monkeypatch):
    monkeypatch.setenv("EXABYTE_USERNAME", "test")
    monkeypatch.setenv("EXABYTE_PASSWORD", "test")
    def login(self):
        return {"X-Account-Id": "test", "X-Auth-Token": "test"}
    monkeypatch.setattr(
        "exabyte_api_client.endpoints.login.LoginEndpoint.login", login
    )
    c = client.Client()
    assert utils.env.exabyte_client_id == "test"
    assert utils.env.exabyte_client_secret == "test"
