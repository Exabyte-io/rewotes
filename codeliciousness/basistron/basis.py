# -*- coding: utf-8 -*-
import bz2
import glob
import os
from collections import defaultdict
from typing import Any, Dict, List, Optional, Tuple

import requests

from basistron import utils

log = utils.get_logger(__name__)


class State:
    """Finite state machine for file parsing."""

    shell_map = {
        "S": 0,
        "P": 1,
        "D": 2,
        "F": 3,
        "G": 4,
        "H": 5,
        "I": 6,
        "J": 7,
        "K": 8,
        "L": 9,
        "M": 10,
    }

    def __init__(self):
        self.live = False
        self.skip = False
        self.scheme = None
        self.symbol = None
        self.shell = 0
        self.func = -1

    def update(self, line: str):
        """Update the state based on contents of a line from
        the NWChem basis set file format."""
        if line.startswith("#BASIS SET:"):
            self.scheme = line.replace("#BASIS SET:", "")
            self.skip = True
        if line.startswith("BASIS"):
            self.live = True
            self.skip = True
        elif line.startswith("END"):
            self.live = False
            self.scheme = None
            self.symbol = None
            self.shell = 0
            self.func = -1

    def emit_tuple(self, line_list: List[str]) -> Optional[Tuple]:
        """Emit structured data when appropriate."""
        try:
            float(line_list[0])
        except ValueError:
            self.symbol = line_list[0]
            self.shell = self.shell_map[line_list[1]]
            self.func += 1
        else:
            return (
                self.scheme,
                self.symbol,
                self.shell,
                self.func,
                float(line_list[0]),
                float(line_list[1]),
            )


class Basis(object):

    BSE_URL = "http://basissetexchange.org/download/current/nwchem/tbz"
    TARBALL = "nwchem_basis_sets.tar.bz2"

    @classmethod
    def download_basis_sets(cls, cache_dir: Optional[str] = None) -> None:
        """This doesn't work quite right so manual download
        instructions are provided in the README.md."""
        cache_dir = cache_dir or utils.default_cache_dir()
        os.makedirs(cache_dir, exist_ok=True)
        resp = requests.get(
            cls.BSE_URL,
            allow_redirects=True,
            headers={"Content-Type": "application/xml"},
        )
        resp.raise_for_status()
        path = os.path.join(cache_dir, cls.TARBALL)
        with bz2.open(path, "w") as f:
            f.write(resp.content)

    @classmethod
    def load_basis_sets(
        cls,
        cache_dir: Optional[str] = None,
        unpacked_dir: Optional[str] = None,
    ) -> Dict[str, Any]:
        """ "Collect all the basis set data for ranking."""
        cache_dir = cache_dir or utils.default_cache_dir()
        unpacked_dir = unpacked_dir or "basis_set_bundle-nwchem-bib"
        basis_dir = os.path.join(cache_dir, unpacked_dir)
        basis_set_files = sorted(glob.glob(os.path.join(basis_dir, "*nw")))
        log.info(f"loading {len(basis_set_files)} basis set files")
        data = {}
        for basis_set_file in basis_set_files:
            try:
                result = cls.parse_basis_set(basis_set_file)
            except KeyError:  # skip support for SP hybrid shells
                continue
            data.update(result)
        return data

    @classmethod
    def rank_basis_sets(cls, basis_data: Dict[str, Any]) -> Dict[str, Any]:
        """Regroups basis sets by symbol and sorts the
        available basis sets by total number of contracted,
        total number of primitive, then number of contracted,
        then number of primitive functions in each shell
        respectively. Assumes basis_data comes from
        load_basis_sets."""
        ordered = defaultdict(list)
        for basis, data in basis_data.items():
            for symbol, contraction in data["contractions"].items():
                ordered[symbol].append(
                    {
                        "basis": basis,
                        "sort_by": cls.get_key_from_contraction(contraction),
                    }
                )
        for sets in ordered.values():
            sets.sort(key=lambda obj: obj["sort_by"])
        return ordered

    @staticmethod
    def get_allowed_basis_sets(ranked: Dict[str, Any], symbols: List[str]) -> List[str]:
        """Gets allowed basis sets for a given set of symbols.
        Maintains sorted order of first symbol provided. Assumes
        ranked takes the form of the output of rank_basis_sets."""
        allowed = None
        for symbol in symbols:
            sets = [obj["basis"] for obj in ranked[symbol]]
            if allowed is None:
                allowed = sets
            allowed = [basis for basis in allowed if basis in sets]
        return allowed

    @staticmethod
    def get_key_from_contraction(
        contraction: str,
    ) -> Tuple[Tuple[int, ...], Tuple[int, ...]]:
        """Create the sort key used for ranking basis sets.
        Updates to the structure for sorting can be extended here."""
        primitive, contracted = contraction.split(" -> ")

        def parse_num(string: str):
            """If it's stupid but it works.."""
            ints, i, s = [], 0, ""
            while i < len(string):
                c = string[i]
                if c.isdigit():
                    s += c
                else:
                    if s:
                        ints.append(int(s))
                    s = ""
                i += 1
            return ints

        return tuple(parse_num(contracted)), tuple(parse_num(primitive))

    @staticmethod
    def parse_basis_set(basis_path: str):
        """Rudimentary nwchem basis set file
        parsing using a finite-state machine."""

        state = State()
        contractions = {}
        data = []

        with open(basis_path, "r") as f:
            for line in f:
                line = line.strip()
                state.update(line)
                if state.live and not state.skip:
                    try:
                        scheme, symbol, *dat = state.emit_tuple(line.split())
                        contractions[symbol] = scheme
                        data.append((symbol, *dat))
                    except TypeError:
                        continue
                state.skip = False

        file_name = basis_path.split("/")[-1]
        basis_name = ".".join(file_name.split(".")[:-1])

        return {
            basis_name: {
                "data": data,
                "contractions": contractions,
            }
        }
