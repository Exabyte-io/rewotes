import os
import importlib
import pytest

from basistron.basis import Basis

class response:
    def __init__(self, text=""):
        self.text = text
        self.status_code = 200
        self.content = b""
    def raise_for_status(self):
        pass

def test_download_basis(monkeypatch, tmppath):
    def get(url, **kwargs):
        return response("basis set data")
    monkeypatch.setattr("basistron.basis.requests.get", get)
    path = (tmppath / "subdir").as_posix()
    Basis.download_basis_sets(cache_dir=path)

def setup_basis(tmppath):
    unpacked_dir = "subdir"
    (tmppath / unpacked_dir).mkdir()
    cache_dir = tmppath.as_posix()
    with open(os.path.join(cache_dir, unpacked_dir, "basis.nw"), "w") as f:
        obj = importlib.resources.read_text("basistron.static", "basis.nw")
        f.write(obj)
    return cache_dir, unpacked_dir

def test_load_basis(tmppath):
    cache_dir, unpacked_dir = setup_basis(tmppath)
    data = Basis.load_basis_sets(cache_dir=cache_dir, unpacked_dir=unpacked_dir)
    assert data

def test_rank_basis(tmppath):
    cache_dir, unpacked_dir = setup_basis(tmppath)
    data = Basis.load_basis_sets(cache_dir=cache_dir, unpacked_dir=unpacked_dir)
    rank = Basis.rank_basis_sets(data)
    assert rank

def test_get_allowed_basis(tmppath):
    cache_dir, unpacked_dir = setup_basis(tmppath)
    data = Basis.load_basis_sets(cache_dir=cache_dir, unpacked_dir=unpacked_dir)
    rank = Basis.rank_basis_sets(data)
    allowed = Basis.get_allowed_basis_sets(rank, ["C", "H", "Kr"])
    assert allowed == ["basis"]
