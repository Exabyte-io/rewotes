# -*- coding: utf-8 -*-

from enum import Enum

from typing import List, Tuple, Optional, Union

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



class Driver(BaseModel):
    """Abstraction layer allowing multiple modes of execution."""
    xyz_data: List[Tuple[str, float, float, float]]
    target_property: str
    reference_value: Optional[Union[str, float]] = None
    reference_tolerance: Optional[float] = 0.01
