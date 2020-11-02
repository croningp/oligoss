import os
import pytest
from oligoss.utils.parameter_handlers.parameter_handlers import generate_parameters
from oligoss.utils.errors import InvalidMSFragmentation

HERE = os.path.abspath(os.path.dirname(__file__))
DATA_FOLDER = os.path.join(HERE, '..', 'data_sets')
INPUTS_FOLDER = os.path.join(HERE, '..', 'input_param_files')

@pytest.mark.unit
def test_instrument_defaults_extractors():
    """
    Tests whether instrument defaults for extractor parameters are filled in
    correctly when instrument-dependent extractor parameters are not explicitly
    stated in input file
    """
    input_json = os.path.join(
        INPUTS_FOLDER,
        "instrument_defaults_extractors.json")

    parameters = generate_parameters(
        params_json=input_json)
    assert "active_instruments" not in parameters.__dict__.keys()
    assert parameters.instrument == "orbitrap_lumos_tribrid"

    assert parameters.extractors.error == 5
    assert parameters.extractors.error_units == "ppm"
    assert parameters.extractors.rt_units == "min"

@pytest.mark.unit
def test_instrument_defaults_silico():
    """
    Tests whether default silico values for parameters dependent on combination
    of instrument AND polymer class are retrieved and correct.
    """
    input_json = os.path.join(
        INPUTS_FOLDER,
        "instrument_polymer_defaults_silico.json")

    parameters = generate_parameters(
        params_json=input_json)
    assert parameters.polymer_class == "depsipeptide"

    assert parameters.silico.ms2.fragment_series == ["b", "y", "a"]
    assert parameters.silico.ms2.min_z == 1
    assert parameters.silico.ms2.max_z == 1
    assert parameters.silico.ms2.max_neutral_losses is None

    assert parameters.silico.ms1.min_z == 1
    assert not parameters.silico.ms1.max_z
    assert parameters.silico.ms1.max_neutral_losses is None

@pytest.mark.unit
def test_invalid_silico_instruemntation():
    """
    Tests whether invalid fragmentation pathways in input files are picked
    up as InvalidMSFragmentation errors
    """
    input_json = os.path.join(
        INPUTS_FOLDER,
        "instrument_polymer_invalid_fragmentation_silico.json")

    with pytest.raises(InvalidMSFragmentation, match=r"valid"):
        parameters = generate_parameters(
            params_json=input_json)
        assert parameters.silico.ms2.fragment_series == ["c", "x"]

@pytest.mark.unit
def test_postprocess_default_instrumentation():

    input_json = os.path.join(
        INPUTS_FOLDER,
        "instrument_polymer_postprocessing_defaults.json")
    parameters = generate_parameters(
        params_json=input_json)

    assert parameters.postprocess.core_linear_series == ["b", "y"]
    assert parameters.postprocess.dominant_signature_cap == 80
    assert parameters.postprocess.optional_core_fragments == ["b1"]
