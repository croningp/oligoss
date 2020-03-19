from silico_utils.monomer import Monomer
from silico_utils.backbone import Backbone
from silico_utils.modification import Modification

from typing import List

class PolymerCompositions:
    def __init__(self, monomers: List[Monomer]):
        self.monomers = monomers

class PolymerFullSilico(PolymerCompositions):
    def __init__(self):
        super().__init__()