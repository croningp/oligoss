from typing import List
from monomer import Monomer

class Polymer(Monomer):

    def __init__(
        self,
        elongation_unit: int,
        min_length: int,
        max_length: int,
        reactivity_classes: List[str, float],
        symmetry: bool,
        **kwargs
    ) -> None:
        super().__init__(locals())
