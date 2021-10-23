# -*- coding: utf-8 -*-

from collections import Counter
from enum import Enum
from typing import Any, Dict, List, Optional, Tuple, Union

from pydantic import BaseModel

from basistron.cccbdb import Cccbdb


class ReferenceRegime(Enum):
    experimental = "experimental"
    calculated = "calculated"

class ExperimentalReferenceProperty(Enum):
    """Map command-line arguments to CCCBDB URIs."""
    polarizability = Cccbdb.EXPT_POL
    vibrational_frequency = Cccbdb.EXPT_VIB
    homo_lumo_gap = Cccbdb.EXPT_IE


class CalculatedReferenceProperty(Enum):
    """Map command-line arguments to CCCBDB URIs."""
    polarizability = Cccbdb.POL_CALC         # units angstrom^3
    vibrational_frequency = Cccbdb.VIB_FREQ  # units cm-1
    homo_lumo_gap = Cccbdb.HOMO_LUMO         # units eV


def validate_property(regime: str, property: str):
    typ = ExperimentalReferenceProperty if (
        regime == ReferenceRegime.experimental.value
    ) else CalculatedReferenceProperty
    try:
        return getattr(typ, property)
    except AttributeError:
        raise Exception(
            f"property {property} not supported in regime {regime}"
        )


class Execution(BaseModel):
    """The state of a given execution."""
    xyz_data: Tuple[Tuple[str, float, float, float], ...]
    property: Union[
        ExperimentalReferenceProperty, CalculatedReferenceProperty,
    ]
    regime: ReferenceRegime
    value: Optional[float] = None
    tolerance: Optional[float] = 1.0

    def acceptable_format(self) -> str:
        if self.value is None:
            return None
        return "{self.value}±{self.tolerance:.2f}%"

    def acceptable_range(self) -> Tuple[float, float]:
        if self.value is None:
            return None
        lower = (1 - self.tolerance) * self.value
        upper = (1 + self.tolerance) * self.value
        return lower, upper

    def simple_formula(self) -> str:
        symbol_count = Counter([r[0] for r in self.xyz_data])
        return ''.join([f"{k}{v}" for k, v in symbol_count.items()])

    def xyz_data_to_dict(self) -> Dict[str, List[Dict[str, Any]]]:
        ang2au = 1.889723
        elements = []
        coordinates = []
        for i, (sym, *val) in enumerate(self.xyz_data):
            i += 1
            elements.append({"id": i, "value": sym})
            coordinates.append(
                {"id": i, "value": [v * ang2au for v in val]}
            )
        return {
            "elements": elements,
            "coordinates": coordinates,
        }
