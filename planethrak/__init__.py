"""PlanetThrak thermal-cracking calculations.

The first development phase reproduces the archived MATLAB lookup-table
workflow before introducing more general fracture constitutive models.
"""

from .fronts import CrackingFront, load_front
from .intersection import CrackingIntersection, find_cracking_intersection

__all__ = [
    "CrackingFront",
    "CrackingIntersection",
    "find_cracking_intersection",
    "load_front",
]
