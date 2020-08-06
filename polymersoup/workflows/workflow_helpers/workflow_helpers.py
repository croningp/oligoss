"""
This file contains functions that call from multiple Polymersoup
modules and / or may be required by several sequencing workflows.
"""

import itertools

from ...extractors.spectra_processing import (
    calculate_ma_consensus_spectra,
    bin_peaks_by_mz,
    find_target_match
)

def calculate_ma_sequence_hits(
    parameters: object,
    ripper_ms2_data: dict,
    confidence_assignments: dict,
    spectral_assignments: dict
) -> dict:

    #  get minimum confidence % for sequences to be assigned MA
    min_confidence = parameters.postprocessing[
        "molecular_assembly"]["min_confidence"]

    #  remove sequences that have not been assigned with sufficient confidence
    #  for calculating MA
    seqs = {
        seq: assignments for seq, assignments in spectral_assignments.items()
        if confidence_assignments[seq] >= min_confidence}
    if not seqs:
        return {}

    #  iterate through sequences and their spectral assignments, working out
    #  how many unique precursors have been assigned to each sequence
    for _ in seqs():

        #  check whether MS2 spectra from multiple unique precursors are to
        #  be combined in consensus spectrum / MA calculation
        if not parameters.postprocessing.molecular_assembly[
                "combine_precursors"]:
            pass

def get_spectra_for_sequence(
    ripper_ms2_data,
    spectral_assignments,
    target_sequence
):

    assignments = spectral_assignments[target_sequence]
    spectrum_ids = spectrum_ids = itertools.chain(
        [list(v) for v in assignments.values()])
    return [
        [spectrum_id, info]
        for spectrum_id, info in ripper_ms2_data.items()
        if spectrum_id in spectrum_ids]


def get_ma_individual_precursors(
    ripper_ms2_data: dict,
    parameters: object,
    min_rel_intensity: float,
    min_abs_intensity: float,
    assignments: dict
) -> dict:

    spectra = get_spectra_for_sequence(
        ripper_ms2_data=ripper_ms2_data,
        spectral_assignments=assignments,
        target_sequence=""
    )

    #  get minimum peak identity
    min_peak_identity = parameters.postprocessing.molecular_assembly[
        "min_peak_identity"]

    #  get list of unique precursors (unique within specified ppm matching
    #   window)
    precursors = bin_peaks_by_mz(
        peaks=[[spectrum["parent"], 1] for spectrum in spectra],
        ppm_window=parameters.postprocessing.ppm_window,
        min_rel_intensity=min_rel_intensity,
        min_abs_intensity=min_abs_intensity,
        intensity_grouping="mean")

    #  iterate through precursors, identifying their MS2 product ion spectra
    for precursor in precursors:
        hit_spectra = [
            x for x in spectra
            if find_target_match(
                target_mass=precursor[0],
                candidate_masses=[x["parent"]],
                error=parameters.postprocessing.ppm_window,
                err_units="ppm")]

        #  generate a separate consensus spectrum for each precursor and
        #  calculate MA of this consensus spectrum
        ma = [calculate_ma_consensus_spectra(
            spectra=hit_spectra,
            ppm_window=parameters.postprocessing,
            intensity_grouping="mean",
            min_abs_intensity=min_abs_intensity,
            min_rel_intensity=min_rel_intensity,
            min_peak_identity=min_peak_identity)]
        return ma
