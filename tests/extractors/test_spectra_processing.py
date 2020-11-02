"""
This file contains tests for spectra processing functions.
"""
import pytest

from oligoss.extractors.spectra_processing import bin_spectrum_peaks

@pytest.fixture
def fake_ms1_spectrum():
    return {
        "100.0000": 10,
        "100.1001": 10000,
        "234.1242": 1200,
        "400.8614": 1400,
        "400.8640": 600,
        "234.1241": 800,
        "234.1251": 200,
        "100.0001": 5000,
        "902.1029": 20,
        "902.1028": 50,
        "mass_list": [
            100.0, 100.0001, 100.1001, 234.1241, 234.1242, 400.8614, 400.8640,
            902.1028, 902.1029
        ],
        "retention_time": 1.61242}


@pytest.mark.unit
def test_ms1_spectrum_binning(fake_ms1_spectrum):
    """
    Tests peak binning function.

    Args:
        fake_ms1_spectrum (dict): spectrum in ripper format.
    """

    #  for random ppm values
    for ppm in [1, 1.1, 3.2, 4.6, 5]:
        five_ppm_bin = bin_spectrum_peaks(
            ppm_window=ppm,
            spectrum=fake_ms1_spectrum,
            min_abs_intensity=50,
            min_rel_intensity=1E-3,
            intensity_grouping="sum"
        )

        test_peaks = [
            [100.1001, 10000],
            [100.0001, 5000],
            [234.1242, 2000],
            [400.8614, 1400],
            [400.864, 600],
            [902.1028, 50]]
        for peak in test_peaks:
            assert five_ppm_bin[str(peak[0])] == peak[1]
