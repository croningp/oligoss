import os
import pytest
from ...utils.parameter_handlers.parameter_handlers import generate_parameters

HERE = os.path.abspath(os.path.dirname(__file__))
DATA_FOLDER = os.path.join(HERE, '..', 'data_sets')
INPUTS_FOLDER = os.path.join(HERE, '..', 'input_param_files')

@pytest.mark.unit
def test_instrument_defaults_orbitrap():

    input_json = os.path.join(
        INPUTS_FOLDER,
        "full_input_parameters_no_defaults.json")

    parameters = generate_parameters(
        params_json=input_json)
    assert "active_instruments" not in parameters.__dict__.keys()
    assert parameters.extractors.error == 5
