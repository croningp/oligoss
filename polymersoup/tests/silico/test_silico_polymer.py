import os
import copy
import pytest

from ...utils.parameter_handlers.parameter_handlers import generate_parameters
from ...utils.errors import InvalidMonomerId
from ...silico.polymer_info.polymer import Polymer, Monomer

HERE = os.path.abspath(os.path.dirname(__file__))
INPUTS_FOLDER = os.path.join(HERE, '..', 'input_param_files')

@pytest.mark.unit
def test_polymer_basic():

    parameters = generate_parameters(
        params_json=os.path.join(
            INPUTS_FOLDER,
            'full_input_parameters_no_defaults.json'))

    fake_polymer_class = copy.deepcopy(parameters)
    fake_polymer_class.polymer_class = "DNA"

    with pytest.raises(Exception, match=r"valid polymer"):
        err_polymer = Polymer(fake_polymer_class)
        for key in err_polymer.__dict__:
            assert key == "params"

    with pytest.raises(InvalidMonomerId):
        polymer = Polymer(parameters)
        for k in polymer.__dict__:
            assert k in ["params", "polymer_config", "monomers"]

    parameters.monomers = ["C", "D", "E"]

    polymer = Polymer(parameters)

    assert not polymer.sidechain_modifications
    assert not polymer.terminal_modifications

    monomer_ids = [monomer for monomer in parameters.monomers]
    polymer_monomer_ids = [
        monomer.id for monomer in polymer.monomers]
    assert monomer_ids == polymer_monomer_ids

    for monomer in polymer.monomers:
        copy_monomer = Monomer(
            monomer_id=monomer.id,
            mode=parameters.mode,
            ms2_signature_types=parameters.silico.ms2.signatures,
            polymer_info=parameters.polymer_class)
        for prop, value in monomer.__dict__.items():
            assert copy_monomer.__dict__[prop] == value
