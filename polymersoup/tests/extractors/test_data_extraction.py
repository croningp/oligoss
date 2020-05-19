import os
import pytest
import mzmlripper.extractor as ripper
from ...extractors.data_extraction import match_mass
from ...extractors.general_functions import open_json

@pytest.fixture
def example_mass_list():
    return {"mass_list": [
            57.1412,
            124.4417,
            722.9747]}

@pytest.mark.unit
def test_mass_match(example_mass_list):
    exact_positive_test = match_mass(
        spectrum=example_mass_list,
        mass_range=[57.1412, 57.1412])

    positive_lower_bound_test = match_mass(
        spectrum=example_mass_list,
        mass_range=[57.1412, 57.1512])

    positive_upper_bound_test = match_mass(
        spectrum=example_mass_list,
        mass_range=[57.1312, 57.1412])

    negative_test = match_mass(
        spectrum=example_mass_list,
        mass_range=[125.4417, 125.5580])

    negative_lower_bound_test = match_mass(
        spectrum=example_mass_list,
        mass_range=[125.4316, 125.4416])

    negative_upper_bound_test = match_mass(
        spectrum=example_mass_list,
        mass_range=[124.4418, 124.4518])

    assert exact_positive_test == ['57.1412']
    assert positive_lower_bound_test == ['57.1412']
    assert positive_upper_bound_test == ['57.1412']
    assert negative_test == []
    assert negative_lower_bound_test == []
    assert negative_upper_bound_test == []

@pytest.mark.unit
def test_mzml_to_json_rt_conversion():
    mzml_folder = os.path.join(
        os.path.dirname(os.path.realpath(__file__)), 'test_mzmls')
    bruker_mzml = os.path.join(mzml_folder, '20190522-EFF MS2 CE5.mzML')
    output_folder = os.path.join(mzml_folder, 'ripper_output')

    not_converted_json = os.path.join(
        output_folder, 'not_converted', 'ripper_20190522-EFF MS2 CE5.json')
    sec_converted_json = os.path.join(
        output_folder, 'sec_converted', 'ripper_20190522-EFF MS2 CE5.json')

    if os.path.isfile(not_converted_json):
        os.remove(not_converted_json)
    if os.path.isfile(sec_converted_json):
        os.remove(sec_converted_json)

    ripper.process_mzml_file(
        filename=bruker_mzml,
        out_dir=os.path.join(output_folder, 'not_converted'))
    not_converted_ripper = open_json(not_converted_json)

    ripper.process_mzml_file(
        filename=bruker_mzml,
        out_dir=os.path.join(output_folder, 'sec_converted'),
        rt_units='sec')
    sec_converted_ripper = open_json(sec_converted_json)

    not_conv_MS1_rts = [
        v['retention_time'] for v in not_converted_ripper['ms1'].values()]
    not_conv_MS2_rts = [
        v['retention_time'] for v in not_converted_ripper['ms2'].values()]

    manual_conv_MS1 = [str(float(rt) / 60) for rt in not_conv_MS1_rts]
    manual_conv_MS2 = [str(float(rt) / 60) for rt in not_conv_MS2_rts]

    sec_conv_MS1_rts = [
        v['retention_time'] for v in sec_converted_ripper['ms1'].values()]
    sec_conv_MS2_rts = [
        v['retention_time'] for v in sec_converted_ripper['ms2'].values()]

    assert manual_conv_MS1.sort() == sec_conv_MS1_rts.sort()
    assert manual_conv_MS2.sort() == sec_conv_MS2_rts.sort()
