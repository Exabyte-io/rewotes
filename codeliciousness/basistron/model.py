# -*- coding: utf-8 -*-

from enum import Enum

from typing import List, Tuple, Optional, Union, Dict, Any

from pydantic import BaseModel


class Property(Enum):
    """Provide different scopes for properties."""

    @classmethod
    def is_valid_property(cls: Enum, prop: str) -> bool:
        """Validate that provided target property is recognized."""
        for sub in cls.__subclasses__():
            if prop in sub.__members__:
                return True
        return False

class SinglePointProperty(Property):
    """Properties only requiring a total energy convergence."""
    energy_convergence = 0
    homo_lumo_gap = 1

class RelaxationProperty(Property):
    """Properties requiring a full relaxation."""
    vibrational_frequencies = 0



class Execution(BaseModel):
    """The state of a given execution."""
    xyz_data: Tuple[Tuple[str, float, float, float], ...]
    target_property: str
    reference_value: Optional[Union[str, float]] = None
    reference_tolerance: Optional[float] = 0.01

    def xyz_data_to_dict(self) -> Dict[str, List[Dict[str, Any]]]:
        elements = []
        coordinates = []
        for i, (sym, *val) in enumerate(self.xyz_data):
            i += 1
            elements.append({"id": i, "value": sym})
            coordinates.append({"id": i, "value": val})
        return {
            "elements": elements,
            "coordinates": coordinates,
        }