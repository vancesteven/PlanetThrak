"""PlanetThrak thermal-cracking calculations.

The first development phase reproduces the archived MATLAB lookup-table
workflow before introducing more general fracture constitutive models.
"""

from .capacity import (
    equiangular_latlon_area_weights,
    fractured_reactive_volume,
    maximum_bound_water_mass,
    spherical_shell_volumes,
)
from .column import FractureColumn, LegacyAccessibility, legacy_front_accessibility
from .field import FractureField
from .fronts import CrackingFront, load_front
from .global_grid import FractureGridResult, fracture_field_from_pt_grids
from .intersection import CrackingIntersection, find_cracking_intersection
from .legacy_structure import RadialColumn, legacy_mars_column
from .planetprofile import fracture_column_from_arrays, fracture_column_from_planetprofile
from .radiogenic import past_radiogenic_heat_uthk

__all__ = [
    "CrackingFront",
    "CrackingIntersection",
    "FractureColumn",
    "FractureField",
    "FractureGridResult",
    "LegacyAccessibility",
    "RadialColumn",
    "equiangular_latlon_area_weights",
    "find_cracking_intersection",
    "fracture_column_from_arrays",
    "fracture_column_from_planetprofile",
    "fracture_field_from_pt_grids",
    "fractured_reactive_volume",
    "legacy_front_accessibility",
    "legacy_mars_column",
    "load_front",
    "maximum_bound_water_mass",
    "past_radiogenic_heat_uthk",
    "spherical_shell_volumes",
]
