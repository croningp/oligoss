from typing import List

class Modification():

    def __init__(
        self,
        mod_mass: float,
        termini: List[int],
        side_chains: List[str],
        free_mod_fragments: dict,
        mod_mass_diff: dict,
        universal_ms2_shift: bool,
        **kwargs
    ) -> None:
        super().__init__(locals())
