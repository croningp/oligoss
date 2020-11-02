import os
import pytest

from oligoss.utils.parameter_handlers.parameter_handlers import generate_parameters

HERE = os.path.abspath(os.path.dirname(__file__))
DATA_FOLDER = os.path.join(HERE, '..', 'data_sets')
INPUTS_FOLDER = os.path.join(HERE, '..', 'input_param_files')


@pytest.mark.unit
def test_full_silico_no_defaults():

    params_json = os.path.join(
        INPUTS_FOLDER,
        "full_input_parameters_no_defaults.json")
    parameters = generate_parameters(
        params_json=params_json,
        param_classes=["core", "silico"])

    assert "core" not in parameters.__dict__.keys()
    assert parameters.screening_method == "exhaustive"
