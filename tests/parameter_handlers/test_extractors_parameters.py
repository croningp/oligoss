import os
import pytest
from oligoss.utils.parameter_handlers.parameter_handlers import generate_parameters
from oligoss.utils.type_dicts.parameter_fallbacks import EXTRACTOR_FALLBACKS

HERE = os.path.abspath(os.path.dirname(__file__))
INPUTS_FOLDER = os.path.join(HERE, '..', 'input_param_files')

@pytest.mark.unit
def test_extractors_basic_typing():
    """
    Tests typing of extractors parameters when wrong types are supplied in
    input file
    """
    input_json = os.path.join(
        INPUTS_FOLDER,
        "full_input_parameters_no_defaults.json")

    parameters = generate_parameters(
        params_json=input_json)

    assert parameters.extractors.error == 5
    assert parameters.extractors.min_ms2_peak_abundance == 100
    assert parameters.extractors.pre_screen_filters[
        "min_ms2_max_intensity"] == 0.3
    assert parameters.extractors.error_units == "ppm"

@pytest.mark.unit
def test_extractors_fallbacks():
    """
    Tests whether default fallback values are used to rescue non-optional
    parameters with missing values
    """
    input_json = os.path.join(
        INPUTS_FOLDER,
        "extractor_fallbacks.json")

    parameters = generate_parameters(
        params_json=input_json)

    assert parameters.extractors.pre_screen_filters

    assert not parameters.extractors.min_ms1_total_intensity
    assert not parameters.extractors.min_ms2_total_intensity
    assert not parameters.extractors.min_ms2_max_intensity

    assert parameters.extractors.min_ms1_max_intensity == 1000

    for param in EXTRACTOR_FALLBACKS:
        assert param in parameters.extractors.__dict__.keys()
