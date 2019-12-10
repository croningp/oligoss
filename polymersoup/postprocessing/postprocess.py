"""
This file contains functions to pass on postprocessing parameters and
execute postprocessing by calling on core functions from postprocess_helpers.py

"""
from .postprocess_helpers import *


def assign_confidence_sequences(
    silico_dict,
    confirmed_dict,
    postprocess_params
):
    """ This function assigns confidence scores to each sequence in the confirmed dict
        and returns a dictionary of sequences and their confidence scores.

    Arguments:
        silico_dict (dict) -- dictionary of all possible sequences (silico dictionary).
        confirmed_dict (dict) -- dictionary of confirmed sequences.
        postprocess_params (dict) -- dictionary of postprocessing parameters, 
            retrieved from input parameters JSON file.

    Returns:
        dict -- dictionary of sequences and their confidence scores.
    """
    # remove sequences from in silico dict that do not have any confirmed
    # fragments
    silico_dict = {
        sequence : subdict
        for sequence, subdict in silico_dict.items()
        if sequence in confirmed_dict
    }

    # initiate dict to store confidence scores for sequences
    confidence_dict = {}

    # iterate through sequences that have confirmed fragments, assigning each
    # a confidence score
    for sequence, confirmed_fragments in confirmed_dict.items():

        print(f"assigning confidence for {sequence}")
        

        # assign sequence a confidence score, using constraints specified in
        # input postprocess parameters from input parameters .json file
        confidence_assignment = assign_confidence_sequence(
            insilico_fragments=silico_dict[sequence]["MS2"],
            confirmed_fragments=confirmed_fragments,
            core_fragment_series=postprocess_params["core_linear_series"],
            optional_fragments=postprocess_params["optional_core_frags"],
            exclude_fragments=postprocess_params["excluded_fragments"],
            essential_fragments=postprocess_params["essential_fragments"],
            essential_fragment_cap=postprocess_params[
                "dominant_signature_cap"],
            sequence_coverage_weight=postprocess_params["subsequence_weight"]
        )

        # add sequence and its final confidence score to confidence_dict
        confidence_dict[sequence] = confidence_assignment

    # return dictionary of sequences and their confidence scores
    return confidence_dict

def find_Rt_I_sequences(
    MS1_EICs,
    MS2_EICs,
    confident_assignments,
    postprocess_parameters,
    confidence_scores,
    return_confidence=True
):
    """
    Takes dicts of MS1 and MS2 EICs, confidently assigned sequences, and fits
    a retention time and intensity to each confidently assigned sequence.

    Args:
        MS1_EICs (dict): dictionary of COMPOSITIONS and corresponding MS1 EICs.
        MS2_EICs (dict): dictionary of SEQUENCES and corresponding MS2 EICs.
        confident_assignments (dict): dictionary of sequences and their
            in silico MS1 and MS2 data for sequences that have been confirmed
            with high confidence.
        postprocess_parameters (dict): dictionary of postprocess parameters
            passed on from input parameters .json file.
        confidence_scores (dict): dictionary of confirmed sequences and 
            assigned confidence scores.
    Kewyword Args:
        return_confidence (bool): specifies whether to return sequences' 
            confidence score in output. Defaults to True.

    Returns:
        dict: dictionary of final intensities and retention times.
    """
    # init dict to store compositions of hit sequences
    compositions = {}

    # iterate through sequences assigned with high confidence
    # and identify their composition in MS1_EICs
    for sequence in confident_assignments:
        sorted_seq = "".join(sorted(sequence))

        composition = [
            seq for seq in MS1_EICs
            if "".join(sorted(seq)) == sorted_seq
        ][0]

        if composition in compositions:
            compositions[composition].append(sequence)
        else:
            compositions[composition] = [sequence]
    final_Rt_Is = {}

    for composition in compositions:

        sequences = [
            seq for seq in compositions[composition]
            if (seq not in MS2_EICs
                and seq in confident_assignments)
        ]

        seq_Rt_Is = get_Rt_I_from_ms1_EIC(
            EIC=MS1_EICs[composition],
            Rt_bin=postprocess_parameters["Rt_bin"],
            backup_Rt_bin=postprocess_parameters["backup_Rt_bin"],
            n_targets=len(sequences),
            min_relative_intensity=postprocess_parameters[
                "min_relative_intensity"]
        )

        print(f'following data points found for {sequences}: {seq_Rt_Is}')
        for i in range(0, len(seq_Rt_Is)):
            final_Rt_Is[sequences[i]] = seq_Rt_Is[i]

        final_Rt_Is.update(
            {
                seq: ['unassigned', 'unassigned']
                for seq in sequences
                if seq not in final_Rt_Is
            }
        )

    # check whether to return sequence confidence scores; if so, add to output
    if return_confidence:
         
        for seq, Rt_I in final_Rt_Is.items():
            print(f"confidence score for {seq} = {confidence_scores[seq]}")
            write_list = [elem for elem in Rt_I]
            write_list.append(confidence_scores[seq])
            print(f"write_list = {write_list}")
            final_Rt_Is[seq] = write_list

    return final_Rt_Is

def postprocess_massdiff_spectral_assignments(
    spectral_assignment_dict,
    ripper_dict,
    ms_level=2
):
    """
    Summarises monomer fingerprint and massdiff assignments to spectra. Summary
    output: N total spectra / N spectra assigned to fingerprint and/or massdiff,
    % assigned spectra in which base peak (dominant ion) has been assigned to
    either a signature ion or found within a mass ladder for one or more 
    target monomers. NOTE: this function is be to used in the 
    mass_difference_screen
    
    Args:
        spectral_assignment_dict (dict): dict of spectrum ids and associated
            assignments for monomer-specific signatures and / or fingerprints.
            For format of this dict, see function: 
            'fingerprint_screen_MSn_from_bpc_precursors' in extractor_helpers.py
        ripper_dict (dict): mass spec data dict in mzml ripper format 
        ms_level (int, optional): specifies MSn level of spectra that have
            been screened and assigned to monomer signatures and / or 
            massdiffs. Defaults to 2.
    
    Returns:
        dict: summary of monomer fingerprint assignments for spectra. Format:
            {
                'N_total
            }
    """

    # filter out relevant msdata from full ripper dict and count number of 
    # recorded MSn spectra (where n = ms_level)
    msdata = ripper_dict[f'ms{ms_level}']
    n_total_spectra = len(msdata.keys())

    # count number of MSn spectra that have been assigned to one or more 
    # monomer-specific fingerprint and / or mass ladder 
    n_assigned_spectra = len([
        key for key in spectral_assignment_dict
        if key.find('spectrum') > -1])

    # count number of assigned MSn spectra with a base peak accounted for by
    # monomer fingerprints and / or mass ladders 
    n_base_peak_assignments = len([
        spectrum for spectrum in spectral_assignment_dict.values()
        if 'base_peak' in spectrum and spectrum['base_peak_assigned']])

    # work out % assigned spectra 
    if n_assigned_spectra > 0: 
        assigned_spectra = (n_assigned_spectra/n_total_spectra)*100
    else:
        assigned_spectra = 0 

    # work out % assigned spectra that have base peak accounted for 
    if n_base_peak_assignments > 0: 
        bp_assignments = (n_base_peak_assignments/n_assigned_spectra)*100 
    else:
        bp_assignments = 0

    # return summary of spectral assignment info 
    return {
        'N_total_spectra': n_total_spectra,
        'total_assigned_spectra': n_assigned_spectra,
        '%_assigned_spectra': assigned_spectra,
        'N_basepeak_assignments': n_base_peak_assignments,
        '%_basepeak_assignments': bp_assignments 
    }
    