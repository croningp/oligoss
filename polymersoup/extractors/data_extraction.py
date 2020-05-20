from .general_functions import os, logging, open_json, write_to_json
from bisect import bisect_left, bisect_right
from .filters import ripper_dict, find_precursors, min_total_intensity_filter

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

        # make sure error is absolute
        if extractor_parameters.error_units == 'ppm':
            error = float(extractor_parameters.error) / 1E6

        else:
            error = float(extractor_parameters.error)

        sequence_EIC = generate_EIC(
            masses=masses,
            ripper_dict=MS1_dict,
            error=error,
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

def match_mass(spectrum, mass_range):
    """ This function takes a mass range for a single target and finds matches
    within the spectrum.

    Args:
        spectrum (dict): single spectrum, must contain 'mass_list'.
        mass_range (list): format: [minimum mass, maximum mass]

    Returns:
        matches (list[float]): list of mass matches from spectrum that fall
        within limit.
    """
    matches = []

    # find closest values to min and max mass in mass list
    mass_list = spectrum['mass_list']
    min_match = bisect_left(mass_list, mass_range[0])
    max_match = bisect_right(mass_list, mass_range[1])

    # find matches within mass range
    try:
        mass_matches = [mass_list[min_match]]

    # if there's no matches (index exceeds the list), return empty list
    except IndexError:
        return matches

    if min_match != max_match:
        mass_matches = mass_list[min_match:max_match]

    if mass_matches:
        matches = [
            format(m, '.4f') for m in mass_matches if (
                m >= mass_range[0]) and (m <= mass_range[1])]

    return matches

def generate_EIC(
        masses, ripper_dict, error, rt_units, min_total_intensity=None):
    """ Takes a list of ions (m/z values) and generates a combined extracted ion
    chromatogram (EIC) of those ions from mass spectra in mzml ripper format
    (ripper_dict).

    Args:
        masses (list[float]): list of target ion m/z values
        ripper_dict (dict): mzml ripper mass spectrometry data dictionary
        error (float): absolute error tolerance for matching target ions to
        masses found in spectra
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

def min_ms2_peak_abundance_filter(
        spectra, peak_list, error, min_ms2_peak_abundance=None):
    """ This function checks all spectra meet the minimum ms2 peak abundance.
    Each spectrum is screened for the most intense peak from the peak list that
    is associated with the sequence. The intensity of this peak is compared with
    the most intense peak in the spectrum (not-related to sequence) and the
    spectrum is returned if it exceeds the min ms2 peak abundance threshold.

    Args:
        spectra (dict): dictionary of MS2 spectra to be filtered.
        peak_list (list(float)): list of all peaks (floats) associated with a
        given sequence.
        error (float): absolute error tolerance for matching target ions to
        masses found in spectra
        min_ms2_peak_abundance (int, optional): minimum percentage of ms2
        peak abundance. Defaults to None.

    Returns:
        min_ms2_peak_filtered (dict): dictionary of spectra that have passed the
        filter.
    """
    # return unfiltered spectra if no filter spectified
    if min_ms2_peak_abundance is None:
        return spectra

    min_ms2_peak_filtered = {}

    # find whether peak associated with sequence passes abundance limit
    for spectrum_id, spectrum_dict in spectra.items():

        # find most intense peak in spectrum
        max_intensity = max(
            [float(x) for x in list(spectrum_dict.values()) if type(x) == int])

        # find most intense peak in spectrum associated to sequence
        highest_int_frag = 0
        for mass in peak_list:
            matches = match_mass(
                spectrum=spectrum_dict,
                mass_range=[mass - error, mass + error])
            for match in matches:
                if spectrum_dict[str(match)] > highest_int_frag:
                    highest_int_frag = spectrum_dict[str(match)]

        ms2_peak_abundance = (
            float(highest_int_frag) / float(max_intensity)) * 100

        # if ms2 peak abundance is higher than threshold, pass the spectrum
        if ms2_peak_abundance >= min_ms2_peak_abundance:
            min_ms2_peak_filtered[spectrum_id] = spectrum_dict

    return min_ms2_peak_filtered

def confirm_fragment(masses, error, spectra):
    """ This function iterates through spectra, looking for mass matches for the
    target masses. A full list of silico masses found is returned.

    Args:
        masses (list[float]): list of target ion m/z values
        error (float): absolute error tolerance for matching target ions to
        masses found in spectra
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

def confirm_all_fragments(fragment_dict, spectra, error):
    """ This function confirms all fragments, including signatures of a silico
    fragment dictionary in a dictionary of MS2 spectra for a single sequence.
    It also creates a json listing which spectra contained fragment hits for
    certain sequences.

    Args:
        fragment_dict (dict): silico MSMS dictionary for the sequence
        spectra (dict): dictionary of MS2 spectra where the precursor mass
            matches the target sequence
        error (float): absolute error tolerance for matching target ions to
            masses found in spectra

    Returns:
        [type]: [description]
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
            spectra=spectra)

        if fragment_matches:
            confirmed_fragments["core"][fragment] = fragment_matches
            all_spectra_matches[fragment] = spectra_matches

    # confirm signatures
    for signature, masses in silico_signatures_list.items():
        signature_matches, sig_spectra_matches = confirm_fragment(
            masses=masses,
            error=error,
            spectra=spectra)

        if signature_matches:
            confirmed_fragments["signatures"][signature] = signature_matches
            all_spectra_matches[signature] = sig_spectra_matches

    # confirm terminal modifications
    for mod, masses in silico_signatures['terminal_modifications'].items():

        mod_matches, mod_spectra_matches = confirm_fragment(
            masses=masses,
            error=error,
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
        error (float): absolute error tolerance for matching target ions to
        masses found in spectra
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
            error=error)

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
                error=error)

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

        # make sure error is absolute
        if extractor_parameters.error_units == 'ppm':
            error = float(extractor_parameters.error) / 1E6
        else:
            error = float(extractor_parameters.error)

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
            error=error,
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
