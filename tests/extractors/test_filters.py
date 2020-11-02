import os
import pytest
import mzmlripper.extractor as ripper
from oligoss.utils.file_io import open_json
from oligoss.extractors.filters import rt_filter, intensity_filter, find_precursor,\
    apply_prefilters_ms2, min_ms2_peak_abundance_filter, match_mass

@pytest.fixture
def example_spectrum():
    return {"spectrum_1": {
            "57.1412": 1743,
            "59.0803": 1819,
            "59.9266": 1895,
            "61.2694": 2052,
            "63.8771": 1910,
            "70.4623": 1871,
            "73.6410": 2196,
            "80.9607": 1827,
            "82.0431": 1980,
            "83.9820": 1687,
            "89.0953": 1934,
            "110.0708": 60556,
            "115.7346": 1779,
            "124.4417": 2020,
            "127.4669": 1957,
            "129.7666": 1829,
            "132.6218": 1985,
            "141.9798": 2230,
            "156.0693": 3181,
            "156.0759": 30415,
            "156.0828": 3590,
            "160.2658": 1999,
            "168.0759": 17091,
            "177.0764": 2027,
            "178.0600": 5909,
            "196.0706": 8613,
            "214.0811": 26594,
            "254.0755": 2739,
            "269.1126": 2706,
            "271.1290": 2164,
            "272.0859": 5600,
            "305.1342": 3904,
            "315.1191": 2387,
            "351.1400": 2959,
            "354.8730": 2107,
            "356.0610": 1946,
            "363.1384": 3349,
            "409.1455": 3947,
            "479.5587": 2272,
            "722.9747": 2146,
            "retention_time": "1.770478261317",
            "mass_list": [
                57.1412,
                59.0803,
                59.9266,
                61.2694,
                63.8771,
                70.4623,
                73.641,
                80.9607,
                82.0431,
                83.982,
                89.0953,
                110.0708,
                115.7346,
                124.4417,
                127.4669,
                129.7666,
                132.6218,
                141.9798,
                156.0693,
                156.0759,
                156.0828,
                160.2658,
                168.0759,
                177.0764,
                178.06,
                196.0706,
                214.0811,
                254.0755,
                269.1126,
                271.129,
                272.0859,
                305.1342,
                315.1191,
                351.14,
                354.873,
                356.061,
                363.1384,
                409.1455,
                479.5587,
                722.9747],
            "parent": "741.2699"}}

@pytest.fixture
def example_mass_list():
    return {"mass_list": [
            57.1412,
            124.4417,
            722.9747]}

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

@pytest.mark.unit
def test_rt_filter(example_spectrum):
    # example spectrum retention time = 1.770478
    positive_test = rt_filter(
        spectra=example_spectrum,
        min_rt=1.75,
        max_rt=1.78)

    negative_higher_test = rt_filter(
        spectra=example_spectrum,
        min_rt=0.75,
        max_rt=1.25)

    negative_lower_test = rt_filter(
        spectra=example_spectrum,
        min_rt=1.90,
        max_rt=2.25)

    no_min_test = rt_filter(
        spectra=example_spectrum,
        min_rt=None,
        max_rt=1.78)

    no_max_test = rt_filter(
        spectra=example_spectrum,
        min_rt=1.75,
        max_rt=None)

    no_filter_test = rt_filter(
        spectra=example_spectrum,
        min_rt=None,
        max_rt=None)

    assert len(positive_test) == 1
    assert len(negative_higher_test) == 0
    assert len(negative_lower_test) == 0
    assert len(no_filter_test) == 1
    assert len(no_min_test) == 1
    assert len(no_max_test) == 1

@pytest.mark.unit
def test_min_max_intensity_filter(example_spectrum):
    # example spectrum max intensity = 60556.0
    positive_test = intensity_filter(
        spectra=example_spectrum,
        min_max_intensity=60000,
        min_total_intensity=None)

    negative_test = intensity_filter(
        spectra=example_spectrum,
        min_max_intensity=60600,
        min_total_intensity=None)

    assert len(positive_test) == 1
    assert len(negative_test) == 0

@pytest.mark.unit
def test_min_total_intensity_filter(example_spectrum):
    # example spectrum total intensity = 230915.0
    positive_test = intensity_filter(
        spectra=example_spectrum,
        min_max_intensity=None,
        min_total_intensity=230000)

    negative_test = intensity_filter(
        spectra=example_spectrum,
        min_max_intensity=None,
        min_total_intensity=250000)

    assert len(positive_test) == 1
    assert len(negative_test) == 0

@pytest.mark.unit
def test_precursor_mass_filter(example_spectrum):
    # example spectrum parent mass = 741.2699
    exact_positive_test_abs = find_precursor(
        spectra=example_spectrum,
        ms2_precursor=741.2699,
        error=0.01,
        error_units='abs')

    exact_positive_test_ppm = find_precursor(
        spectra=example_spectrum,
        ms2_precursor=741.2699,
        error=5,
        error_units='ppm')

    positive_test = find_precursor(
        spectra=example_spectrum,
        ms2_precursor=741.2730,
        error=0.01,
        error_units='abs')

    negative_test_lower_bound = find_precursor(
        spectra=example_spectrum,
        ms2_precursor=741.2589,
        error=0.01,
        error_units='abs')

    negative_test_higher_bound = find_precursor(
        spectra=example_spectrum,
        ms2_precursor=741.2800,
        error=0.01,
        error_units='abs')

    assert len(exact_positive_test_abs) == 1
    assert len(exact_positive_test_ppm) == 1
    assert len(positive_test) == 1

    assert len(negative_test_lower_bound) == 0
    assert len(negative_test_higher_bound) == 0

@pytest.mark.unit
def test_min_ms2_peak_abundance_filter(example_spectrum):
    # example spectrum max intensity = 60556.0
    exact_positive_test = min_ms2_peak_abundance_filter(
        spectra=example_spectrum,
        peak_list=[110.0708],
        error=0.01,
        min_ms2_peak_abundance=100)

    positive_test = min_ms2_peak_abundance_filter(
        spectra=example_spectrum,
        peak_list=[156.0759],
        error=0.01,
        min_ms2_peak_abundance=50)

    no_filter_test = min_ms2_peak_abundance_filter(
        spectra=example_spectrum,
        peak_list=[156.0759],
        error=0.01,
        min_ms2_peak_abundance=None)

    negative_test = min_ms2_peak_abundance_filter(
        spectra=example_spectrum,
        peak_list=[141.9798],
        error=0.01,
        min_ms2_peak_abundance=100)

    assert len(exact_positive_test) == 1
    assert len(positive_test) == 1
    assert len(no_filter_test) == 1
    assert len(negative_test) == 0

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
def test_apply_prefilters_ms2(example_spectrum):

    positive_test = apply_prefilters_ms2(
        spectra=example_spectrum,
        min_rt=1,
        max_rt=2,
        min_max_intensity=60000,
        min_total_intensity=230000,
        ms2_precursors=["741.2699"],
        error=0.01,
        error_units='abs')

    full_negative_test = apply_prefilters_ms2(
        spectra=example_spectrum,
        min_rt=2,
        max_rt=3,
        min_max_intensity=80000,
        min_total_intensity=250000,
        ms2_precursors=["743.2699"],
        error=0.01,
        error_units='abs')

    assert len(positive_test) == 1
    assert len(full_negative_test) == 0
