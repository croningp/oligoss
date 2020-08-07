import logging

from functools import lru_cache

from .filters import (
    find_precursors,
    min_ms2_peak_abundance_filter,
    match_mass
)


def extract_MS1_EICs(
    MS1_silico,
    extractor_parameters,
    MS1_dict,
    filename
):
    """ This function extracts EICs for each sequence in the MS1 silico dict
    from the MS1_ripper_dict.

    Args:
        MS1_silico (Dict[str, List[float]]): full MS1 silico dictionary of the
            format: {'sequence' :[masses]}
        extractor_parameters (object): full extractor parameters
        MS1_dict (Dict[str, dict]): MS1 ripper mass spectrometry data
        filename (str): str name of sample (ripper json name without .json)

    Returns:
        Dict[str, List[List[float]]]: dict of compositions / sequences that have
            met the required MS1 ion intensity threshold. Format:
            {"sequence": [[rt, I]...] ...} where rt, I = retention time and I =
            sequence intensity at rt.
    """

    # format sequence sequence: [[rt, I], ..]
    retention_times = retrieve_retention_times(ripper_dict=MS1_dict)
    ms1_EICs = {}

    if len(MS1_silico) > 1:

        #  log number of compositions being screened for
        logging.info(f"generating MS1 EICs for {len(MS1_silico)} compositions")

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
            ms1_EICs[sequence] = sequence_EIC
            logging.info(f"{sequence} MS1 hit found")

    #  add list of retention times for all ms1 spectra for saving to output
    ms1_EICs['retention_times'] = retention_times

    return ms1_EICs

def retrieve_retention_times(ripper_dict):
    """ This function generates a list of all retention times from a ripperdict.
    This list will be useful in plotting EICs by providing full time series.

    Args:
        ripper_dict (dict): ripper mass spectrometry data dictionary

    Returns:
        retention times (List[float]): list of all retention times
    """
    retention_times = []

    for spectrum in ripper_dict.values():
        retention_times.append(spectrum['retention_time'])

    return retention_times

def generate_EIC(
    masses,
    ripper_dict,
    error,
    error_units,
    rt_units,
    min_total_intensity=None
):
    """ Takes a list of ions (m/z values) and generates a combined extracted ion
    chromatogram (EIC) of those ions from mass spectra in mzml ripper format
    (ripper_dict).

    Args:
        masses (List[float]): list of target ion m/z values
        ripper_dict (Dict[str, dict]): mzml ripper mass spectrometry data
            dictionary
        error (float): error tolerance for matching target ions to
            masses found in spectra
        error_units (str): 'ppm' or 'abs'
        rt_units (str): units of the mzml ripper retention times (min or sec)
        min_total_intensity (float, optional): maximum total intensity of EIC
            for it to be returned. Defaults to None

    Returns:
        EIC (List[List[float]]: extracted ion chromatogram in format:
            [[rt, I]..] where rt = retention time (float) and I = intensity of
            all species that match target ions at that retention time.
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
        masses (List[float]): list of target ion m/z values
        error (float): tolerance for matching target ions to
            masses found in spectra
        error_units (str): 'ppm' or 'abs'
        spectra (Dict[str, Dict[str, float]]): dictionary of all MS2 spectra for
            a single ripper

    Returns:
        matches (List[float]): list of masses found in spectra that match
            original masses
        spectra_matches (Dict[str, List(float)]): dict of spectrum id's where
            fragments were found and which masses were matched
    """
    matches = []
    spectra_matches = {}

    for spectrum_id, spectrum in spectra.items():
        error_value = error
        for mass in masses:

            # make sure error is in absolute units
            if error_units == 'ppm':
                error_value = (float(mass) / 1E6) * error

            # search for matches in spectrum
            mass_range = [mass - error_value, mass + error_value]

            match = match_mass(
                spectrum=spectrum,
                mass_range=mass_range)

            if match:
                matches.append(mass)
                spectra_matches[spectrum_id] = match

    # sort matches by ascending number
    matches = sorted(list(set(matches)))

    return matches, spectra_matches

def confirm_ms2_fragments_isomeric_sequences_concurrent(
    composition,
    silico_dict,
    ripper_ms2,
    params,
    ripper_tag
):

    #  get spectra with matching precursors for target composition
    candidate_spectra = find_precursors(
        ms2_precursors=silico_dict[composition]["MS1"],
        spectra=ripper_ms2,
        error=params.extractors.error,
        error_units=params.extractors.error_units
    )
    if not candidate_spectra:
        return None, None

    candidate_spectra = min_ms2_peak_abundance_filter(
        spectra=candidate_spectra,
        peak_list=silico_dict[composition]["MS2"]["peak_list"],
        error=params.extractors.error,
        min_ms2_peak_abundance=params.extractors.min_ms2_peak_abundance
    )
    if not candidate_spectra:
        return None, None

def confirm_all_fragments_concurrent(
    fragment_dict,
    ms2_spectra,
    ripper_tag,
    params
):
    """
    Takes MS/MS silico dict for isomeric sequences and screens MS2 spectra for
    target sequences.

    Args:
        fragment_dict (Dict[str, dict]): dict of sequences and corresponding
            MS1 precursor and MS2 fragment info.
        ms2_spectra (List[dict]): list of ripper spectra.
        ripper_tag (str): unique str tag for original ripper file. NOTE: this
            is only used in caching, but is essential when screening multiple
            rippers in multiprocessing workflows.
        params (Parameters): parameters object.

    Returns:
        dict, dict: confirmed fragments, spectral matches
    """
    # separate out 'core' fragments from signatures
    silico_fragments = {
        k: v for k, v in fragment_dict["MS2"].items()
        if k not in ["composition", "signatures"]
    }
    silico_signatures = fragment_dict["MS2"]["signatures"]
    silico_signatures_list = {
        k: v for k, v in silico_signatures.items()
        if k != 'terminal_modifications'}
    confirmed_fragments = {
        "core": {}, "signatures": {"terminal_modifications": {}}}
    all_spectra_matches = {}

    #  get spectra with precursor matches
    ms2_spectra = find_precursors(
        ms2_precursors=fragment_dict["MS1"],
        spectra=ms2_spectra,
        error=params.extractors.error,
        error_units=params.extractors.error_units
    )
    if not ms2_spectra:
        return None, None

    #  make sure spectra meet MAPI criteria
    ms2_spectra = min_ms2_peak_abundance_filter(
        spectra=ms2_spectra,
        error=params.extractors.error,
        peak_list=fragment_dict["peak_list"],
        min_ms2_peak_abundance=params.extractors.min_ms2_peak_abundance
    )

    @lru_cache(maxsize=1000)
    def confirm_masses_internal(masses, spectrum_ids, tag):
        """
        Inner function for searching MS2 spectra for fragment masses. NOTE:
        this is a memoized function. tag MUST be assigned to the correct
        source RipperDict.

        Args:
            masses (Tuple[float]): target m/z values.
            spectrum_ids (Tuple[str]): spectrum_ids for spectra matching
                precursors and MAPI. NOTE: used only for caching.
            tag (str): unique string tag for source ripper. NOTE: used purely
                for caching.

        Returns:
            List[float]: list of fragments with match in target spectra.
        """
        return confirm_fragment(
            masses=masses,
            error=params.extractors.error,
            error_units=params.extractors.error_units,
            spectra=ms2_spectra
        )

    #  get spectrum ids for caching
    spectrum_ids = tuple(sorted([spectrum_id for spectrum_id in ms2_spectra]))

    # confirm 'core' fragments
    for fragment, masses in silico_fragments.items():
        fragment_matches, spectra_matches = confirm_masses_internal(
            masses=tuple(masses),
            spectrum_ids=spectrum_ids,
            tag=ripper_tag)

        if fragment_matches:
            confirmed_fragments["core"][fragment] = fragment_matches
            all_spectra_matches[fragment] = spectra_matches

    # confirm signatures
    for signature, masses in silico_signatures_list.items():
        signature_matches, sig_spectra_matches = confirm_masses_internal(
            masses=tuple(masses),
            spectrum_ids=spectrum_ids,
            tag=ripper_tag)

        if signature_matches:
            confirmed_fragments["signatures"][signature] = signature_matches
            all_spectra_matches[signature] = sig_spectra_matches

    # confirm terminal modifications
    if "terminal_modifications" in silico_signatures:
        for mod, masses in silico_signatures["terminal_modifications"].items():

            mod_matches, mod_spectra_matches = confirm_masses_internal(
                masses=tuple(masses),
                spectrum_ids=spectrum_ids,
                tag=ripper_tag)

            if mod_matches:
                confirmed_fragments[
                    "signatures"]["terminal_modifications"][mod] = mod_matches
                all_spectra_matches[mod] = mod_spectra_matches

    return confirmed_fragments, all_spectra_matches

def confirm_ms2_fragments_isomeric_sequences(
    target_composition,
    precursors,
    full_msms_silico_dict,
    ripper_data,
    params
):
    """
    Takes a target composition string, a full MS/MS silico dict, ripper data,
    and then screens for all potential sequences that match composition of
    target.

    Args:
        target_composition (str): composition string for target isomers
        precursors (List[float]): list of precursor m/z values for composition
        full_msms_silico_dict (Dict[str, dict]): dict of sequences and
            associated MS1 precursors, MS2 product ion subdicts
        ripper_data (RipperDict): RipperDict object
        params (Parameters): Parameters object

    Returns:
        Dict[str, dict], Dict[str, dict]: confirmed_fragments, spectral_matches
    """

    #  log target composition
    logging.info(
        "searching for MS2 fragments for sequences isomeric to\n"
        f"{target_composition}")
    confirmed_fragment_dict, confirmed_spectral_matches = {}, {}

    candidate_spectra = find_precursors(
        ms2_precursors=precursors,
        spectra=ripper_data.ms2,
        error=params.extractors.error,
        error_units=params.extractors.error_units
    )

    if not candidate_spectra:
        return confirmed_fragment_dict, confirmed_spectral_matches

    sequence_entries = [
        seq for seq in [x for x in full_msms_silico_dict if x != "compositions"]
        if (full_msms_silico_dict[seq]["MS2"]["composition"]
            == target_composition)
    ]

    #  log number of isomers to be screened at MS2
    logging.info(
        f"screening {len(sequence_entries)} for composition\n"
        f"{target_composition}. Searching {len(candidate_spectra)} MS2 spectra")

    for sequence in sequence_entries:

        # make sure that minimum ms2 peak abundance is reached
        filtered_precursor_spectra = min_ms2_peak_abundance_filter(
            spectra=candidate_spectra,
            peak_list=full_msms_silico_dict[sequence]["peak_list"],
            error=params.extractors.error,
            min_ms2_peak_abundance=(
                params.extractors.min_ms2_peak_abundance
            )
        )

        # search for MS2 fragment hits in filtered spectra
        if filtered_precursor_spectra:

            confirmed_fragments, spectra_matches = confirm_all_fragments(
                fragment_dict=full_msms_silico_dict[sequence]["MS2"],
                spectra=filtered_precursor_spectra,
                error=params.extractors.error,
                error_units=params.extractors.error_units
            )
            confirmed_fragment_dict[sequence] = confirmed_fragments
            confirmed_spectral_matches[sequence] = spectra_matches

    #  log number of sequences isomeric to target have been found at MS2
    logging.info(
        f"found {len(confirmed_fragment_dict)} sequence isomeric to\n"
        f"{target_composition} with one or more confirmed MS2 fragments")

    return confirmed_fragment_dict, confirmed_spectral_matches

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
        confirmed_fragments (Dict[str, List[float]]): dictionary of all
            confirmed fragments where key is the fragment string and values
            are a list of silico masses
        all_spectra_matches (Dict[str, List[str]]): dictionary of all spectra in
            which fragments were confirmed. format:
            {'fragment string': ['spectra id's (str)']}
    """
    # separate out 'core' fragments from signatures
    silico_fragments = {
        k: v for k, v in fragment_dict.items()
        if k not in ["composition", "signatures"]
    }
    silico_signatures = fragment_dict['signatures']
    silico_signatures_list = {
        k: v for k, v in silico_signatures.items()
        if k != 'terminal_modifications'}
    confirmed_fragments = {
        "core": {}, "signatures": {"terminal_modifications": {}}}
    all_spectra_matches = {}

    @lru_cache(maxsize=1000)
    def confirm_fragment_internal(masses):
        """
        This inner function is a memoized version of the confirm_fragment
        function. This is to save repeatedly screening for isomeric fragments
        when screening isomeric sequences. Should save time. NOTE: please be
        careful with cache size and reduce if it takes up too much memory.
        Args:
            masses (Tuple): tuple of fragment m/z values.
        """

        return confirm_fragment(
            masses=masses,
            error=error,
            error_units=error_units,
            spectra=spectra
        )

    # confirm 'core' fragments
    for fragment, masses in silico_fragments.items():
        fragment_matches, spectra_matches = confirm_fragment_internal(
            tuple(masses))

        if fragment_matches:
            confirmed_fragments["core"][fragment] = fragment_matches
            all_spectra_matches[fragment] = spectra_matches

    # confirm signatures
    for signature, masses in silico_signatures_list.items():
        signature_matches, sig_spectra_matches = confirm_fragment_internal(
            tuple(masses))

        if signature_matches:
            confirmed_fragments["signatures"][signature] = signature_matches
            all_spectra_matches[signature] = sig_spectra_matches

    # confirm terminal modifications
    if "terminal_modifications" in silico_signatures:
        for mod, masses in silico_signatures["terminal_modifications"].items():

            mod_matches, mod_spectra_matches = confirm_fragment_internal(
                tuple(masses))

            if mod_matches:
                confirmed_fragments[
                    "signatures"]["terminal_modifications"][mod] = mod_matches
                all_spectra_matches[mod] = mod_spectra_matches

    #  clear cache for inner function. NOTE: not clearing the cache WILL result
    #  in false results if calling this function to screen more than one ripper
    confirm_fragment_internal.cache_clear()

    return confirmed_fragments, all_spectra_matches
