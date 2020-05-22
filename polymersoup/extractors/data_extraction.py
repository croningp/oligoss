from .general_functions import os, logging, open_json, write_to_json
from .filters import ripper_dict, find_precursors, min_total_intensity_filter,\
    min_ms2_peak_abundance_filter, match_mass

logging.basicConfig(
    format='%(message)s - %(asctime)s',
    datefmt='%H:%M:%S %m/%d/%Y ',
    level=logging.INFO)

def extract_MS1_EICs(
        MS1_silico,
        extractor_parameters,
        MS1_dict,
        filename,
        output_folder):
    """ This function extracts EICs for each sequence in the MS1 silico dict
    from the MS1_ripper_dict.

    Args:
        MS1_silico (dict): full MS1 silico dictionary of the format:
        {'sequence' :[masses]}
        extractor_parameters (object): full extractor parameters
        MS1_dict (dict): MS1 ripper mass spectrometry data
        filename (str): str name of sample (ripper json name without .json)
        output_folder (str): full filepath to output folder

    Returns:
        trimmed_silico (dict): MS1 silico dictionary whereby any sequence not
        found in sufficient abundance to build an EIC, is removed
        MS1_fragments (dict): format: {'sequence': [[rt, I]]}
    """
    # format sequence sequence: [[rt, I], ..]
    retention_times = retrieve_retention_times(ripper_dict=MS1_dict)
    MS1_EIC = {}
    MS1_fragments = {}
    trimmed_silico = MS1_silico.copy()

    # generate EIC for each
    for sequence, masses in MS1_silico.items():

        sequence_EIC = generate_EIC(
            masses=masses,
            ripper_dict=MS1_dict,
            error=extractor_parameters.error,
            error_units=extractor_parameters.error_units,
            min_total_intensity=extractor_parameters.min_ms1_total_intensity,
            rt_units=extractor_parameters.rt_units)

        if sequence_EIC:
            logging.info(f'EIC generated for {sequence}')
            MS1_EIC[sequence] = sequence_EIC

            # sort composition alphabetically
            MS1_fragments[(''.join(sorted(sequence)))] = masses

        # if no EIC was generated for the sequence, remove it from silico dict
        if not sequence_EIC:
            logging.info(f'{sequence} not in sufficient abundance for EIC')
            del trimmed_silico[sequence]

    MS1_EIC['retention_times'] = retention_times

    # write MS1 EIC to json
    write_to_json(MS1_EIC, os.path.join(
        output_folder, f'{filename}_MS1_EICs.json'))
    return trimmed_silico, MS1_fragments

def retrieve_retention_times(ripper_dict):
    """ This function generates a list of all retention times from a ripperdict.
    This list will be useful in plotting EICs by providing full time series.

    Args:
        ripper_dict (dict): ripper mass spectrometry data dictionary

    Returns:
        retention times (list[float]): list of all retention times
    """
    retention_times = []

    for spectrum in ripper_dict.values():
        retention_times.append(spectrum['retention_time'])

    return retention_times

def generate_EIC(
        masses, ripper_dict, error, error_units, rt_units,
        min_total_intensity=None):
    """ Takes a list of ions (m/z values) and generates a combined extracted ion
    chromatogram (EIC) of those ions from mass spectra in mzml ripper format
    (ripper_dict).

    Args:
        masses (list[float]): list of target ion m/z values
        ripper_dict (dict): mzml ripper mass spectrometry data dictionary
        error (float): error tolerance for matching target ions to
        masses found in spectra
        error_units (str): 'ppm' or 'abs'
        rt_units (str): units of the mzml ripper retention times (min or sec)
        min_total_intensity (float, optional): maximum total intensity of EIC
        for it to be returned. Defaults to None

    Returns:
        EIC (list(list)): extracted ion chromatogram in format: [[rt, I]..]
        where rt = retention time (float) and I = intensity of all species that
        match target ions at that retention time.
    """
    EIC = []
    total_intensity = 0

    for mass in masses:

        # make sure error is in absolute units
        if error_units == 'ppm':
            error = (float(mass) / 1E6 * error)

        # define minimum and maximum mass match considering errors
        mass_range = [mass - error, mass + error]

        for spectrum in ripper_dict.values():

            # if there are matches within range, add rt and I to EIC
            matches = match_mass(spectrum=spectrum, mass_range=mass_range)

            if matches:
                for match in matches:
                    EIC.append(
                        [spectrum['retention_time'],
                            spectrum[str(match)]])
                    total_intensity += spectrum[str(match)]

    # check total intensity meets minimum value, if not return empty EIC
    if min_total_intensity:
        if (total_intensity <= min_total_intensity):
            EIC = []

    return EIC

def confirm_fragment(masses, error, error_units, spectra):
    """ This function iterates through spectra, looking for mass matches for the
    target masses. A full list of silico masses found is returned.

    Args:
        masses (list[float]): list of target ion m/z values
        error (float): tolerance for matching target ions to
        masses found in spectra
        error_units (str): 'ppm' or 'abs'
        spectra (dict): dictionary of all MS2 spectra for a single ripper

    Returns:
        matches (list[float]): list of masses found in spectra that match
        original masses
        spectra_matches (list[str]): list of spectrum id's where fragments were
        found
    """
    matches = []
    spectra_matches = []

    for spectrum_id, spectrum in spectra.items():
        for mass in masses:

            # make sure error is in absolute units
            if error_units == 'ppm':
                error = (float(mass) / 1E6 * error)

            # search for matches in spectrum
            mass_range = [mass - error, mass + error]
            match = match_mass(
                spectrum=spectrum,
                mass_range=mass_range)
            if match:
                matches.append(mass)
                spectra_matches.append(spectrum_id)

    # sort matches by ascending number
    matches = list(set(matches))
    matches.sort()
    spectra_matches = list(set(spectra_matches))
    return matches, spectra_matches

def confirm_all_fragments(fragment_dict, spectra, error, error_units):
    """ This function confirms all fragments, including signatures of a silico
    fragment dictionary in a dictionary of MS2 spectra for a single sequence.
    It also creates a json listing which spectra contained fragment hits for
    certain sequences.

    Args:
        fragment_dict (dict): silico MSMS dictionary for the sequence
        spectra (dict): dictionary of MS2 spectra where the precursor mass
            matches the target sequence
        error (float): error tolerance for matching target ions to
            masses found in spectra
        error_units (str): 'ppm' or 'abs'

    Returns:
        confirmed_fragments (dict): dictionary of all confirmed fragments where
        key is the fragment string and values are a list of silico masses
        all_spectra_matches (dict): dictionary of all spectra in which fragments
        were confirmed. format: {'fragment string': ['spectra id's (str)']}
    """
    # separate out 'core' fragments from signatures
    silico_fragments = {
        k: v for k, v in fragment_dict['MS2'].items() if k != 'signatures'}
    silico_signatures = fragment_dict['MS2']['signatures']
    silico_signatures_list = {
        k: v for k, v in silico_signatures.items()
        if k != 'terminal_modifications'}
    confirmed_fragments = {
        "core": {}, "signatures": {"terminal_modifications": {}}}
    all_spectra_matches = {}

    # confirm 'core' fragments
    for fragment, masses in silico_fragments.items():
        fragment_matches, spectra_matches = confirm_fragment(
            masses=masses,
            error=error,
            error_units=error_units,
            spectra=spectra)

        if fragment_matches:
            confirmed_fragments["core"][fragment] = fragment_matches
            all_spectra_matches[fragment] = spectra_matches

    # confirm signatures
    for signature, masses in silico_signatures_list.items():
        signature_matches, sig_spectra_matches = confirm_fragment(
            masses=masses,
            error=error,
            error_units=error_units,
            spectra=spectra)

        if signature_matches:
            confirmed_fragments["signatures"][signature] = signature_matches
            all_spectra_matches[signature] = sig_spectra_matches

    # confirm terminal modifications
    for mod, masses in silico_signatures['terminal_modifications'].items():

        mod_matches, mod_spectra_matches = confirm_fragment(
            masses=masses,
            error=error,
            error_units=error_units,
            spectra=spectra)

        if mod_matches:
            confirmed_fragments[
                "signatures"]["terminal_modifications"][mod] = mod_matches
            all_spectra_matches[mod] = mod_spectra_matches

    return confirmed_fragments, all_spectra_matches

def confirm_all_sequences_fragments(
        MSMS_insilico_dict,
        MS2_spectra,
        error,
        error_units,
        filename,
        output_folder,
        min_ms2_total_intensity,
        min_ms2_peak_abundance):
    """ This function interates through each sequence in the MSMS_insilico
    dictionary and searches for matches in spectra that mass the sequence
    precursor mass. A dictionary containing sequences, their confirmed fragments
    and the spectra in which these fragments were found is also created. This
    will be used by the postprocessing module to plots spectral assignments.

    Args:
        MSMS_insilico_dict (dict): silico MSMS dictionary for all
        sequences
        MS2_spectra (dict): all MS2 spectra for the ripper file
        error (float): error tolerance for matching target ions to
        masses found in spectra
        error_units (str): 'abs' or 'ppm'
        filename (str): str name of sample (ripper json name without .json)
        output_folder (str): filepath to where confirmed fragment dict will be
        saved as a json
        min_ms2_total_intensity (int): minimum intensity for the total
        intensity for each MS2 spectrum
        min_ms2_peak_abundance (int): minimum intensity percentage for intensity
        of a peak associated with the sequences compared to the most intense
        peak in the spectrum

    Returns:
        confirmed_fragment_dict (dict): dictionary of all sequences with
        confirmed fragments format = {sequence: [confirmed fragment strings]}
    """
    confirmed_fragment_dict = {}
    all_spectra_matches = {}

    # iterate through sequences in MSMS insilico dictionary
    for sequence, frag_dict in MSMS_insilico_dict.items():

        # find MS2 spectra that match precursor
        precursor_spectra = find_precursors(
            spectra=MS2_spectra,
            ms2_precursors=frag_dict['MS1'],
            error=error,
            error_units=error_units)

        # make sure spectra meet min MS2 total intensity threshold
        int_filtered_precursor_spectra = min_total_intensity_filter(
            spectra=precursor_spectra,
            min_total_intensity=min_ms2_total_intensity)

        # make sure that minimum ms2 peak abundance is reached
        filtered_precursor_spectra = min_ms2_peak_abundance_filter(
            spectra=int_filtered_precursor_spectra,
            peak_list=frag_dict['peak_list'],
            error=error,
            min_ms2_peak_abundance=min_ms2_peak_abundance)

        # search for MS2 fragment hits in filtered spectra
        if filtered_precursor_spectra:

            confirmed_fragments, spectra_matches = confirm_all_fragments(
                fragment_dict=frag_dict,
                spectra=filtered_precursor_spectra,
                error=error,
                error_units=error_units)

            if confirmed_fragments:
                confirmed_fragment_dict[sequence] = confirmed_fragments

                all_spectra_matches[sequence] = spectra_matches

    # save confirmed fragment dict and all spectra matches to jsons
    write_to_json(confirmed_fragment_dict, os.path.join(
        output_folder, f'{filename}_confirmed_MS2_fragment_dict.json'))
    write_to_json(all_spectra_matches, os.path.join(
        output_folder, f'{filename}_MS2_spectra_matches.json'))

    return confirmed_fragment_dict

def standard_extraction(MS1_silico, rippers, extractor_parameters, output):
    """ This function iterates through each ripper file, confirming MS1
    compositions from the MS1_silico dict, MS1 EICs are generated.
    Full MSMS insilico fragments are generated for confirmed compositions.
    MS2 fragments are confirmed for each sequence.

    Args:
        MS1_silico (dict): theoretical dictionary of all compositions (k) and
            a list of associated their masses (v)
        rippers (list): list of full filepaths to ripper jsons
        extractor_parameters (object): full extractor parameters
        output (str): full filepath to desired output folder

    Returns:
        extracted_data_folders (list[str]): full filepath to extracted data
        folder for all rippers.
    """
    extracted_data_folders = []
    for ripper_file in rippers:

        # open ripper json
        ripper_data = open_json(ripper_file)
        ripper = ripper_dict(ripper_data)

        # create output folder
        ripper_name = (ripper_file.split('\\')[-1]).replace('.json', '')
        ripper_output = os.path.join(output, ripper_name)
        extracted_data_folders.append(ripper_output)
        if not os.path.exists(ripper_output):
            os.mkdir(ripper_output)

        logging.info(f'extracting MS1 EICs for {ripper_name}')

        # extract MS1 EICs
        trimmed_MS1_silico, MS1_fragments = extract_MS1_EICs(
            MS1_silico=open_json(MS1_silico),
            extractor_parameters=extractor_parameters,
            MS1_dict=ripper.ms1,
            filename=ripper_name,
            output_folder=ripper_output)

        # generate MSMS insilico for confirmed MS1 compositions and save to json
        # UPDATE TO DAVID'S FUNCTION, use trimmed_MS1_silico
        MSMS_insilico = {}
        logging.info(
            'MSMS insilico file location needed to be added to data extraction')
        logging.info(f'confirming MS2 fragments for {ripper_name}')

        # confirm fragments for all sequences
        confirmed_fragment_dict = confirm_all_sequences_fragments(
            MSMS_insilico_dict=MSMS_insilico,
            MS2_spectra=ripper.ms2,
            error=extractor_parameters.errorerror,
            error_units=extractor_parameters.error_units,
            filename=ripper_name,
            output_folder=ripper_output,
            min_ms2_peak_abundance=extractor_parameters.min_ms2_peak_abundance,
            min_ms2_total_intensity=extractor_parameters.min_ms2_total_intensity
        )

        # save full dict of confirmed sequence fragments (MS1 & MS2) to json
        all_confirmed_fragments = {}

        for sequence, ms2_frags in confirmed_fragment_dict.items():

            # sort sequence alphabetically to find confirmed MS1 matches
            sorted_seq = ''.join(sorted(sequence))
            ms1_frags = MS1_fragments[sorted_seq]
            all_confirmed_fragments[sequence] = {
                'MS1': ms1_frags,
                'MS2': ms2_frags}
        confirmed_fragments_output = os.path.join(
            ripper_output, f'{ripper_name}_confirmed_fragment_dict.json')
        write_to_json(all_confirmed_fragments, confirmed_fragments_output)

    return extracted_data_folders
