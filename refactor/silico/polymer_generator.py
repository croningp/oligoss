from backbone import Backbone
from monomer import Monomer

class Polymer(Backbone, Monomer):

    def __init__(self):
        super().__init__(locals())
