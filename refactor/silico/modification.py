from typing import List

class Modification():

    def __init__(
        self,
        modification_id: str,
        mod_mass: float,
        termini: List[int],
        side_chains: List[str],
        free_mod_fragments: dict,
        mod_mass_diff: dict,
        universal_ms2_shift: bool
    ) -> None:
        pass
