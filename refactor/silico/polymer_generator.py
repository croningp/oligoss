from silico_utils.monomer import Monomer

from typing import List

class PolymerCompositions:
    def __init__(self, monomers: List[Monomer]):
        self.monomers = monomers

class PolymerFullSilico(PolymerCompositions):
    def __init__(self):
        super().__init__()
