import logging

from functools import lru_cache

from itertools import chain

from .filters import (
    find_precursors,
    min_ms2_peak_abundance_filter,
    match_mass
)

from ..utils.run_utils import exception_handler
from ..postprocessing.postprocess_helpers import list_unconfirmed_fragments

@exception_handler(verbose=True)
def extract_MS1_EICs(
    MS1_silico,
    extractor_parameters,
    MS1_dict,
    filename,
    db
):
    """
    This function extracts EICs for each sequence in the MS1 silico dict
    from the MS1_ripper_dict.

    Args:
        MS1_silico (Dict[str, List[float]]): full MS1 silico dictionary of the
            format: {'sequence' :[masses]}
        extractor_parameters (object): full extractor parameters
        MS1_dict (Dict[str, dict]): MS1 ripper mass spectrometry data
        filename (str): str name of sample (ripper json name without .json)
        db (bool): specifies whether EICs are being generated from standard
            ripper objects or ripper data pulled from MongoDB database. Defaults
            to False.

    Returns:
        Dict[str, List[List[float]]]: dict of compositions / sequences that have
            met the required MS1 ion intensity threshold. Format:
            {"sequence": [[rt, I]...] ...} where rt, I = retention time and I =
            sequence intensity at rt.
    """
    ms1_EICs = {}

    if len(MS1_silico) > 1:

        #  log number of compositions being screened for
        logging.info(f"generating MS1 EICs for {len(MS1_silico)} compositions")

    # generate EIC for each
    for sequence, masses in MS1_silico.items():

        #  check whether data to be screened is from ripper object or retrieved
        #  from MongoDB database and extract accordingly
        if db:
            sequence_EIC = generate_EIC_db(
                masses=masses,
                ripper_dict=MS1_dict,
                error=extractor_parameters.error,
                error_units=extractor_parameters.error_units,
                min_total_intensity=(
                    extractor_parameters.min_ms1_total_intensity),
                rt_units=extractor_parameters.rt_units
            )
        else:
            sequence_EIC = generate_EIC(
                masses=masses,
                ripper_dict=MS1_dict,
                error=extractor_parameters.error,
                error_units=extractor_parameters.error_units,
                min_total_intensity=(
                    extractor_parameters.min_ms1_total_intensity),
                rt_units=extractor_parameters.rt_units)

        if sequence_EIC:
            ms1_EICs[sequence] = sequence_EIC
            logging.info(f"{sequence} MS1 hit found")
        else:
            logging.info(f"{sequence} not found at MS1")

    return ms1_EICs

def extract_MS1_EICs_db(
    MS1_silico,
    extractor_parameters,
    ripper_name,
    filename,
    connection
):
    """
    This function extracts EICs for each sequence in the MS1 silico dict
    from the MS1_ripper_dict.

    Args:
        MS1_silico (Dict[str, List[float]]): full MS1 silico dictionary of the
            format: {'sequence' :[masses]}
        extractor_parameters (object): full extractor parameters
        ripper_name (str): ripper_name (str): unique tag for ripper file.
            Used for caching and locating correct database connections.
        filename (str): str name of sample (ripper json name without .json)
        connection (str): mongodb connection string.

    Returns:
        Dict[str, List[List[float]]]: dict of compositions / sequences that have
            met the required MS1 ion intensity threshold. Format:
            {"sequence": [[rt, I]...] ...} where rt, I = retention time and I =
            sequence intensity at rt.
        maximum intensity (float): maximum intensity value for sequence.
    """

    # retrieve ms1 ripper data (spectra only, not _id field)
    MS1_dict = list(
        connection.polymersoup[
            f'{ripper_name}_ms1_data'].find({}, {"_id": 0, "spectrum": 1}))

    ms1_EICs = {}

    #  log number of compositions being screened for
    if len(MS1_silico) > 1:
        logging.info(f"generating MS1 EICs for {len(MS1_silico)} compositions")

    # generate EIC for each
    max_intensity = 0

    for sequence, masses in MS1_silico.items():

        sequence_EIC = generate_EIC_db(
            masses=masses,
            ripper_dict=MS1_dict,
            error=extractor_parameters.error,
            error_units=extractor_parameters.error_units,
            min_total_intensity=extractor_parameters.min_ms1_total_intensity,
            rt_units=extractor_parameters.rt_units)

        if sequence_EIC:

            logging.info(f"{sequence} MS1 hit found")
            ms1_EICs[sequence] = sequence_EIC
            max_seq_int = max([seq[1] for seq in sequence_EIC])

            if max_seq_int > max_intensity:
                max_intensity = max_seq_int

    #  log summary of MS1 extraction
    logging.info(
        f"EICs generated for {len(ms1_EICs) - 1} compositions")

    return ms1_EICs, max_intensity

def retrieve_retention_times(ripper_dict):
    """ This function generates a list of all retention times from a ripperdict.
    This list will be useful in plotting EICs by providing full time series.

    Args:
        ripper_dict (dict): ripper mass spectrometry data dictionary

    Returns:
        retention times (List[float]): list of all retention times
    """
    retention_times = []

    for spectrum in ripper_dict:
        retention_times.append(ripper_dict[spectrum]['retention_time'])
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

def generate_EIC_db(
    masses,
    ripper_dict,
    error,
    error_units,
    rt_units,
    min_total_intensity=None
):
    """ Takes a list of ions (m/z values) and generates a combined extracted ion
    chromatogram (EIC) of those ions from mass spectra in mzml ripper format
    (ripper_dict). NOTE: this function is only for generating EICs from ripper
    data pulled from MongoDB database.

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

        for spec in ripper_dict:

            spectrum = spec["spectrum"]

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

def confirm_all_fragments_concurrent(
    fragment_dict,
    ms2_spectra,
    params
):
    """
    Takes MS/MS silico dict for isomeric sequences and screens MS2 spectra for
    target sequences. NOTE: to be used in multiprocessed workflows.

    Args:
        fragment_dict (Dict[str, List[float]]): dict of fragments and
            corresponding m/z values.
        ms2_spectra (List[Dict[str, dict]]): list of ripper spectra.
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
    def confirm_masses_internal(masses, spectrum_id):
        """
        Inner function for searching MS2 spectra for fragment masses. NOTE:
        this is a memoized function. When screening fragment masses and
        spectra, ensure masses and spectrum_ids are passed in as tuples (or
        other immutable, hashable data structure) and in a consistent order.

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
            spectrum_id=spectrum_ids)

        if fragment_matches:
            confirmed_fragments["core"][fragment] = fragment_matches
            all_spectra_matches[fragment] = spectra_matches

    # confirm signatures
    for signature, masses in silico_signatures_list.items():
        signature_matches, sig_spectra_matches = confirm_masses_internal(
            masses=tuple(masses),
            spectrum_id=spectrum_ids)

        if signature_matches:
            confirmed_fragments["signatures"][signature] = signature_matches
            all_spectra_matches[signature] = sig_spectra_matches

    return confirmed_fragments, all_spectra_matches


def confirm_all_fragments_concurrent_db(
    composition,
    precursors,
    params,
    connection,
    ripper_name
):
    """
    Takes MS/MS silico dict for isomeric sequences and screens MS2 spectra for
    target sequences. NOTE: this function uses MongoDB to access ripper data
    and store screening results.

    Args:
        composition (str): composition string.
        precursors (List[float]): list of precursor masses.
        params (Parameters): parameters object.
        connection (str): polymersoup localhost connection string
        ripper_name (str): unique str tag for original ripper file. NOTE: this
            is only used in caching, but is essential when screening multiple
            rippers in multiprocessing workflows.
    """
    # retrieve ms2 silico information for composition
    polymersoupdb = connection["polymersoup"]

    # retrieve all entries with a compositional match as a list
    ms2_silico = list(polymersoupdb[f"{ripper_name}_ms2_silico"].find(
        {"composition": composition}))

    # set up mongodb collection to store ms2 hit information
    ms2_hits = polymersoupdb[f"{ripper_name}_ms2_hits"]

    # close connection once ms2 silico dict is retrieved and collections made
    connection.close()

    # log target composition
    logging.info("searching for MS2 fragments for sequences isomeric to\n"
                 f"{composition}")

    #  get spectra from ripper ms2_data collection with precursor matches
    ms2_spectra = find_precursors(
        ms2_precursors=precursors,
        error=params.extractors.error,
        error_units=params.extractors.error_units,
        connection=connection,
        ripper_name=ripper_name)

    if not ms2_spectra:
        return None, None

    # log number of isomers to be screened at MS2
    logging.info(
        f"screening for {len(ms2_silico)} isomers for composition\n"
        f"{composition}. Searching {len(ms2_spectra)} MS2 spectra")

    #  iterate through isomeric seqs, getting any confirmed fragments,
    # unconfirmed fragments and their spectral matches
    confirmed_isomer_count = 0

    for dic in ms2_silico:
        sequence = dic["_id"]

        silico_dict = {k: v for k, v in dic.items() if k in [
            "signatures", "series"]}

        # make list of all MS2 masses associated with the sequence
        frag_list = list(
            dic["series"].values()) + list(dic["signatures"].values())

        # un-nest values to make one big list of all peaks associated with that
        # sequence
        peak_list = list(chain.from_iterable(frag_list))

        #  make sure spectra meet MAPI criteria
        ms2_filtered = min_ms2_peak_abundance_filter(
            spectra=ms2_spectra,
            error=params.extractors.error,
            peak_list=peak_list,
            min_ms2_peak_abundance=params.extractors.min_ms2_peak_abundance)

        # search for ms2 fragment hits in filtered spectra
        if ms2_filtered:

            confirmed_frags, spectral_matches = confirm_all_fragments(
                fragment_dict=silico_dict,
                spectra=ms2_filtered,
                error=params.extractors.error,
                error_units=params.extractors.error_units)

            unconfirmed_fragments = list_unconfirmed_fragments(
                confirmed_fragments=confirmed_frags,
                silico_dict=silico_dict)

            # add confirmed fragments and spectral matches to collections
            # tagged by composition
            if confirmed_frags:
                confirmed_isomer_count += 1
                ms2_hits.insert(
                    {
                        "_id": sequence,
                        "confirmed_fragments": confirmed_frags,
                        "unconfirmed_fragments": unconfirmed_fragments,
                        "composition": composition,
                        "spectral_matches": spectral_matches
                    }
                )

    # log number of sequences isomeric to target that have been found at MS2
    logging.info(
        f"found {confirmed_isomer_count} sequences isomeric to \n"
        f"{composition} with one or more confirmed MS2 fragments")

    connection.close()

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
    silico_fragments = fragment_dict["series"]
    silico_signatures = fragment_dict['signatures']
    silico_signatures_list = {
        k: v for k, v in silico_signatures.items()
        if k != 'terminal_modifications'}
    confirmed_fragments = {
        "core": {}, "signatures": {"terminal_modifications": {}}}
    all_spectra_matches = {}

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

    return confirmed_fragments, all_spectra_matches
