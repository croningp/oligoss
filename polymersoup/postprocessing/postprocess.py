import os
import logging
import pandas as pd
from ..extractors.general_functions import open_json, write_to_json
from ..silico.helpers.helpers import get_composition
from .postprocess_helpers import (
    assign_confidence_sequences,
    assign_isomeric_groups,
    all_spectral_assignment_plots
)

def postprocess_ripper(ripper_folder, postprocess_parameters, ms2_data):
    """ This function performs all standard postprocessing operations for a
        single ripper output.

    Args:
        ripper_folder (str): filepath to extracted data folder for a single
        ripper(where all output files from data extraction were saved).
        postprocess_parameters (obj): postprocessing parameters.
        ms2_data (dict): MS2 data only from ripper.

    """

    all_ssw = postprocess_parameters.subsequence_weight
    ripper_name = ripper_folder.split('\\')[-1]

    # retrieve full silico MSMS dict from extracted data
    full_silico_MSMS = open_json(
        os.path.join(
            ripper_folder,
            [file for file in os.listdir(ripper_folder)
                if 'PRE_fragment_screening_silico_dict' in file][0]))

    # retrieve confirmed fragment dict from extracted data
    confirmed_fragments = open_json(
        os.path.join(
            ripper_folder,
            [file for file in os.listdir(ripper_folder)
                if 'confirmed_fragment_dict' in file][0]))

    # retrieve MS1 EICs from extracted data
    MS1_EICs = open_json(
        os.path.join(
            ripper_folder,
            [file for file in os.listdir(ripper_folder)
                if 'MS1_EICs' in file][0]))

    # retrieve spectra matches from extracted data
    MS2_spectra_matches = open_json(
        os.path.join(
            ripper_folder,
            [file for file in os.listdir(ripper_folder)
                if 'spectra_matches' in file][0]))

    # get confidence scores for sequences at all subsequence weights
    for ssw in all_ssw:
        confidence_scores = assign_confidence_sequences(
            silico_dict=full_silico_MSMS,
            confirmed_dict=confirmed_fragments,
            postprocess_params=postprocess_parameters,
            ssw=ssw)

        # create new folder directory for current subsequence weight
        confidence_dir = os.path.join(
            ripper_folder,
            'postprocessing',
            f'confidence_assignments_{ssw}ssw')

        # spectral assignment plots over confidence threshold
        if postprocess_parameters.spectral_assignment_plots:
            plot_sequence_list = [
                s for s, c in confidence_scores.items()
                if float(c) >= postprocess_parameters.min_plot_confidence]

            all_spectral_assignment_plots(
                sequence_list=plot_sequence_list,
                MS_data=ms2_data,
                MS2_spectra_matches=MS2_spectra_matches,
                postprocess_folder=confidence_dir)

        if not os.path.exists(confidence_dir):
            os.makedirs(confidence_dir)

        # write confidence scores to json
        write_to_json(write_dict=confidence_scores, output_json=os.path.join(
            confidence_dir, f'confidence_scores_{ssw}_ssw.json'))

        # convert confidence scores dict into dataframe with columns
        # 'Sequence' and 'Confidence'
        postprocess_summary = pd.DataFrame(
            list(confidence_scores.items()),
            columns=['Sequence', 'Confidence'])
        sequences = list(postprocess_summary['Sequence'])

        # create 'Confirmed Core Fragments' column
        postprocess_summary['Confirmed Core Fragments'] = [
            list(confirmed_fragments[sequence]["core"].keys())
            for sequence in sequences]

        # create 'Confirmed Signatures' column
        postprocess_summary['Confirmed Signatures'] = [
            list(confirmed_fragments[sequence]["signatures"].keys())[:-1]
            for sequence in sequences]

        # create 'Max Intensity' column
        postprocess_summary['Max Intensity'] = [
            max(MS1_EICs[get_composition(sequence)], key=lambda x:x[1])[1]
            for sequence in sequences]

        # create 'Isomeric_Group' column by matching MS1 fragments
        isomeric_groups = assign_isomeric_groups(
            sequences=sequences)
        postprocess_summary['Isomeric_Group'] = postprocess_summary[
            'Sequence'].map(isomeric_groups)

        # write postprocess summary to csv
        postprocess_summary.to_csv(
            os.path.join(confidence_dir, f'postprocess_summary_{ssw}ssw.csv'),
            index=False)

    logging.info(f'postprocessing complete for {ripper_name}')
