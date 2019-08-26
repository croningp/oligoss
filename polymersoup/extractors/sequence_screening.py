"""
This file contains functions for extracting sequence-specific data from
spectra

"""
from .run_extractors import *

def generate_EIC(
    ions,
    ms_level,
    ripper_dict,
    err,
    err_abs=True,
    min_max_intensity=None,
    min_total_intensity=None,
    precursors=None,
    mapi=0,
    peak_list=None
):
    """
    Takes a list of ions (m/z values) and generates a combined extracted ion
    chromatogram (EIC) of those ions from mass spectra in mzml ripper format
    (ripper_dict)

    Args:
        ions (list of floats): list of target ion m/z values
        ms_level (int): MS level at which to screen for ions
        ripper_dict (dict): mzml ripper mass spectrometry data dictionary
        err (float): error tolerance for matching target ions to masses found
                in spectra; units either in absolute mass units or ppm
        err_abs (bool, optional): specifies whether err is in absolute mass
                units or ppm; if True, err is in absolute mass units.
                Defaults to True.
        min_max_intensity (float, optional): minimum maximum intensity of EIC
                for it to be returned. Defaults to None.
        min_total_intensity (float, optional): maximum total intensity of EIC
                for it to be returned. Defaults to None.
        precursors (list of floats, optional): list of precursor ions for
                matching to parents/precursors of MSn spectra
        mapi (optional, float): minimum annotated peak intensity (in % of most
                intense peak), specifies minimum relative intensity of most
                intense peak associated with a sequence for spectra to be used
                in screening
        peak_list (optional, list of floats): list of ALL ions associated with
                targets, not just ions used in generating EIC
    Returns:
        EIC (list of lists): extracted ion chromatogram in format: [[Rt, I]..]
                where Rt = retention time (float) and I = intensity of all
                species that match target ions at that retention time
    """

    # extract spectra at target ms level
    ripper_dict = ripper_dict[f'ms{ms_level}']

    # initiate list to store EIC
    EIC = []

    # iterate through spectra, and list masses found in each spectrum
    for spectrum_info in ripper_dict.values():

        # convert masses to floats
        masses = [float(mass) for mass in spectrum_info['mass_list']]

        # assume spectrum is suitable for searching targets
        valid_spectrum = True

        # if precursors (i.e. 'parents') are specified, check whether precursor
        # / parent of this spectrum matches one or more of the precursors
        if precursors:
            parent = float(spectrum_info['parent'])
            parent_matches = find_multiple_targets(
                targets=precursors,
                candidates=[parent],
                err=err,
                err_abs=err_abs
            )
            # if the spectrum parent does not match proposed precursors, the
            # spectrum is not suitable for screening
            if not parent_matches:
                valid_spectrum = False

        # check whether a minimum annotated peak intensity has been specified;
        # if so, check whether spectrum exceeds this minimum annotated peak
        # intensity
        if mapi > 0:

            # check for a peak_list; if none has been input, mapi cannot be
            # checked
            if not peak_list:
                raise Exception("mapi check cannot be carried out if peak_list is not supplied")

            # calculate annotated peak intensity (in % of most intense peak)
            if mapi > find_maximum_annotated_peak_intensity(
                spectrum=spectrum_info,
                err=err,
                err_abs=err_abs,
                peak_list=peak_list
            ):
                # if annotated peak intensity < mapi, spectrum is not suitable
                # for screening
                valid_spectrum = False

        # initiate list for storing masses that match target ions
        matches = []

        # check whether spectrum is suitable before looking for matches
        if valid_spectrum:

            # iterate through targets, and add matches as they are found for each
            # target to matches list
            for ion in ions:
                matches.extend(find_target(
                    target=ion,
                    candidates=masses,
                    err=err,
                    err_abs=err_abs)
                )

        # get intensity of all matches in spectrum
        intensity = sum_intensity_targets_spectrum(
            spectrum_dict=spectrum_info,
            found_targets=matches
        )
        if matches and intensity == 0:
            raise Exception(f'no intensity found for {matches}')
        # retrieve spectrum's retention time
        retention_time = float(spectrum_info['retention_time'])

        # add retention time, intensity as a sublist to EIC
        if intensity > 0:
            EIC.append([retention_time, intensity])

    # sort EIC by intensity
    # Rt_I = sublists of [retention_time, intensity] in EIC
    EIC.sort(key = lambda Rt_I: Rt_I[1])

    # check for specified minimum max intensity threshold. If this is specified
    # and most intense point on EIC does not exceed this, return blank EIC
    if min_max_intensity:
        if EIC:
            if EIC[-1][1] < min_max_intensity:
                EIC = []

    # check for specified minimum total intensity threshold. If this is
    # specified and sum of EIC intensity does not exceed this, return blank EIC
    if min_total_intensity:
        if EIC:
            if sum([Rt_I[1] for Rt_I in EIC]) < min_total_intensity:
                EIC = []

    # return EIC
    return EIC

def generate_MS1_EICs_sequence_dict(
    silico_dict,
    ripper_dict,
    err,
    err_abs,
    min_max_intensity,
    min_total_intensity
):
    """
    Takes a silico_dict of sequences and corresponding MS1 ions, and outputs
    MS1 EICs for sequences

    Args:
        silico_dict ([type]): [description]
        ripper_dict ([type]): [description]
        err ([type]): [description]
        err_abs ([type]): [description]
        min_max_intensity ([type]): [description]
        min_total_intensity ([type]): [description]
    """
    for sequence in silico_dict:
        if type(silico_dict[sequence]) == dict:
            masses = silico_dict[sequence]["MS1"]
        else:
            masses = silico_dict[sequence]

        silico_dict[sequence] = [float(mass) for mass in masses]

        sequence_MS1_EIC = generate_EIC(
            masses,
            1,
            ripper_dict,
            err,
            err_abs,
            min_max_intensity,
            min_total_intensity
        )

def confirm_fragments_sequence(
    precursors,
    fragment_dict,
    peak_list,
    ripper_dict,
    err,
    err_abs=True,
    ms_level=2,
    min_annotated_peak_assignment=90,
    essential_signatures=[]
):
    """
    This function takes a list of precursor ions, fragment dictionary and
    returns list of confirmed fragments associated with precursors

    Args:
        precursors (list of floats): list of m/z values for sequence precursors
        fragment_dict (dict): dictionary of fragments and associated masses
            in insilico fragment dict format
        peak_list (list of floats): list of ALL MS1 and MS2 / MSn masses
            associated with sequence
        ripper_dict (dict): dict of mass spec data in mzml ripper format
        err (float): error tolerance used when matching masses; units can
            either be in absolute mass units or ppm
        err_abs (bool, optional): specifies whether err units are in absolute
            mass units or ppm; if True, absolute mass unit is used.
            Defaults to True.
        ms_level (int, optional): specifies MS level for screening fragments
            and precursors. Defaults to 2.
        min_annotated_peak_assignment (int, optional): minimum % intensity of
            most intense peak in a given spectrum that matches a mass in the
            peak list. Defaults to 90.
        essential_signatures (list, optional): list of signatures that MUST
            be found in a spectrum for fragments to be confirmed from that
            spectrum. Defaults to [].
    Returns:
        confirmed_fragments (list of strings): list of fragment ids for
            fragments in fragment_dict that have been confirmed as present in
            mass spectrum data from ripper_dict
    """
    # remove spectra at fragmentation ms level that do not have matching
    # precursors
    filtered_dict = filter_parent_ion(
        msdata=ripper_dict,
        parent_ions=precursors,
        err=err,
        err_abs=err_abs,
        ms_level=ms_level)
    
    # pull out spectra at fragmentation ms level
    filtered_dict = filtered_dict[f'ms{ms_level}']

    # initiate list to store confirmed fragment ids
    confirmed_fragments = []

    # iterate through spectra
    for spectrum_info in filtered_dict.values():

        # calculate annotated_peak_assignemnt from spectrum, peak_list
        annotated_peak_assignment = find_maximum_annotated_peak_intensity(
            spectrum=spectrum_info,
            peak_list=peak_list,
            err=err,
            err_abs=err_abs
        )

        # if annotated_peak_assignment is below specified lower threshold,
        # pass over spectrum without searching for fragments
        if annotated_peak_assignment < min_annotated_peak_assignment:
            pass

        # if annotated_peak_assignment is equal to or exceeds threshold, search
        # spectrum for fragments provided in fragment_dict
        else:

            # add fragments found in spectrum to list of confirmed fragments
            # associated with sequence
            confirmed_fragments.extend(
                find_fragments_spectrum(
                    spectrum=spectrum_info,
                    fragment_dict=fragment_dict,
                    err=err,
                    err_abs=err_abs,
                    essential_signatures=essential_signatures
                )
            )
    # remove any duplicates from confirmed fragment list
    confirmed_fragments = list(set(confirmed_fragments))

    # finally, return list of confirmed fragments
    return confirmed_fragments

def find_most_intense_peak_spectrum(
    spectrum
):
    """
    Takes a spectrum in mzml_ripper format and returns tuple of most intense
    peak with m/z and associated intensity in format (m/z, I) where m/z =
    m/z value of most intense peak and I = intensity of most intense peak

    Args:
        spectrum (dict): single spectrum dictionary in mzml ripper format

    Returns:
        most_intense (tuple): tuple of peak m/z and intensity in format
                    (m/z, I)
    """


    masses = []

    # get list of mass (m/z), intensity tuples for all peaks in spectrum
    for mass in spectrum['mass_list']:
        mass_key = str(mass) + "0000"
        mass_key = mass_key[0:mass_key.find('.')+5]
        masses.append(
            (float(mass), float(spectrum[mass_key]))
        )

    # sort peaks by most intense
    masses = sorted(masses, key = lambda x: x[1])

    if type(masses) != list:
        raise Exception(f'lambda fucked')

    if type(masses) != list:
        raise Exception(f'lambda fucked')
    # retrieve most intense m/z in format (m/z, I) where m/z = peak m/z and I
    # = intensity of peak with m/z
    most_intense = masses[-1]

    return most_intense

def find_most_intense_matching_peak(
    spectrum,
    peak_list,
    err,
    err_abs):
    """
    takes a list of theoretical peaks associated with a sequence, a spectrum
    containing observed peaks and associated intensities, and returns the
    intensity of the most intense matching peak

    Args:
        spectrum (dict): spectrum dict in mzml ripper dict format
        peak_list (list): list of theoretical m/z values associated with a
                        sequence
        err (float):
        err_abs ([type]): [description]

    Returns:
        top_intensity (float): intensity of most intense matching peak
    """

    # set top matching intensity to 0; this will be reset if matching peaks
    # are found with intensity > 0
    top_intensity = 0
  
    intensities = []

    # iterate through theoretical peaks associated with sequence
    for peak_mass in peak_list:

        # find real peaks in spectrum that match theoretical sequence peak
        # mass
        matches = find_target(
            target=peak_mass,
            candidates=[float(mass) for mass in spectrum['mass_list']],
            err=err,
            err_abs=err_abs
        )
    
        # get intensities of all matching peaks, sorted in ascending order
        for match in matches:
            try:
                intensities.append(float(spectrum[f'{match}']))
            except KeyError:
                match = str(match) + "0000"
                match = match[0:match.find(".")+5]
                intensities.append(float(spectrum[f'{match}']))
    
    intensities = sorted(intensities)
    
    if not intensities:
        return top_intensity
    else:
        if intensities[-1] > top_intensity:
            top_intensity = intensities[-1]

    return top_intensity

def find_maximum_annotated_peak_intensity(
    spectrum,
    peak_list,
    err,
    err_abs
):
    """
    Searches a spectrum in mzml ripper dict format for peaks (m/z values) in
    a peak list and returns % intensity of most intense matching peak relative
    to most intense peak in the spectrum

    Args:
        spectrum (dict): spectrum dict in mzml ripper format
        peak_list (list of floats): list of peaks (m/z values), typically peaks
                    associated with a target sequence
        err (float): error tolerance in matching peak masses - can be in
                    absolute mass units or ppm
        err_abs (bool): specifies whether err units are in absolute mass units
                    or ppm; if True, units = absolute mass units

    Returns:
        float: maximum annotated peak intensity in % of most intense peak
    """

    # find intensity of most intense peak in spectrum
    most_intense_peak = find_most_intense_peak_spectrum(spectrum)[-1]

    # find intensity of most intense spectrum peak that matches something in
    # the peak list
    most_intense_matching_peak = find_most_intense_matching_peak(
        spectrum=spectrum,
        peak_list=peak_list,
        err=err,
        err_abs=err_abs
    )
   
    # return maximum annotated peak intensity in % of most intense peak in
    # spectrum
    if most_intense_matching_peak > 0:
        return (most_intense_matching_peak/most_intense_peak) * 100
    else:
        return 0

def find_fragments_spectrum(
    spectrum,
    fragment_dict,
    err,
    err_abs,
    essential_signatures=[]
):
    """
    Takes an MS2 / MSn fragment dict, MS2 / MSn spectrum and returns list of
    fragments that have matching ions in the spectrum

    Args:
        spectrum (dict): spectrum dict in mzml ripper format
        fragment_dict (dict): fragment dictionary in standard insilico format
        err (float): error tolerance for matching fragment masses; units can
                    either be absolute mass units or ppm
        err_abs (bool): specifies whether err units are in absolute mass
                    units or ppm; if True, err is in absolute mass units
        essential_signtures (list): list of signature fragment ids that must
                    be found in spectrum if any fragments are to be confirmed
                    from said spectrum (defaults to [])

    Returns:
        confirmed_fragments (list of strings): list of fragment ids from
                    input fragment_dict that have one or more matching masses
                    in the spectrum
    """
    # initiate list to store confirmed fragments
    confirmed_fragments = []

    # retrieve masses from spectrum
    spectrum_masses = [float(mass) for mass in spectrum['mass_list']]

    # iterate through SIGNATURE FRAGMENTS
    for signature, signature_masses in fragment_dict['signatures'].items():
        
        # retrieve masses in spectrum that match fragment target masses
        matches = find_multiple_targets(
            signature_masses,
            spectrum_masses,
            err,
            err_abs
        )

        # if any spectrum masses are a match for fragment masses, add
        # fragment id to list of confirmed fragments
        if matches:
            confirmed_fragments.append(signature)

    # check whether there are any essential signatures, i.e. signatures that
    # MUST be found in the spectrum for fragments to be confirmed
    if essential_signatures:

        # if any essential signatures are not found, return empty list of
        # confirmed fragments
        for signature in essential_signatures:
            if signature not in confirmed_fragments:
                return []

    # iterate through NON-SIGNATURE fragments
    for fragment, masses in fragment_dict.items():

        matches = []
        if fragment != 'signatures':
            # retrieve masses in spectrum that match fragment target masses
            matches = find_multiple_targets(
                masses,
                spectrum_masses,
                err,
                err_abs
            )

        # if any spectrum masses are a match for fragment masses, add
        # fragment id to list of confirmed fragments
        if matches:
            confirmed_fragments.append(fragment)

    # return list of fragments found in spectrum
    return confirmed_fragments

def confirm_fragments_sequence_dict(
    silico_dict,
    ripper_dict,
    extractor_parameters
):
    """
    Takes a full MSMS in silico sequence dict, mass spec data in mzml ripper
    dict format, and returns a dictionary of sequences and MS2 fragments that
    can be confirmed as associated with those sequences from the ripper data

    Args:
        silico_dict (dict): full MSMS in silico sequence dict in format:
            {
                seq: {
                    "MS1": [m/z...],
                    "MS2": {
                        "frag": [m/z...],
                        "signatures": {
                            "signature": [m/z...]
                    },
                    "peak_list": [m/z...]
                    }
                }
            }
        where seq = sequence string; m/z = ion m/z values; "frag" = fragment id
        string; "signature" = fragment id string for monomer-specific
        signature fragment; "peak_list" = list of ALL MS1 AND MS2 ions
        associated with sequence, including signature ions

        ripper_dict (dict): mass spec data dict in mzml ripper format
        err (float): error tolerance for screening masses
        err_abs (bool): specifies whether units of err are in absolute mass
        units or ppm; if True, units = absolute mass units
        min_ms2_peak_abundance (float): minimum relative abundance (in % of
        most abundant / intense MS2 spectrum peak) of most intense peak in
        peak_list
        essential_signatures (list of strings): list of signatures that
        MUST be found in a spectrum for fragments to be confirmed from that
        spectrum
    """

    # check whether extractor_parameters is a string; if so, it will be the
    # file path to input parameters .json file; open extractor parameters as
    # dict
    if type(extractor_parameters) == str:
        extractor_parameters = open_json(
            extractor_parameters
            )["extractor_parameters"]

    # inititate dict for storing sequences and their confirmed fragments, if
    # any are confirmed
    confirmed_fragment_dict = {}

    # iterate through sequences in the in silico MSMS dict
    for sequence, subdict in silico_dict.items():

        # retrieve confirmed fragments for sequence from ripper_dict
        confirmed_fragments = confirm_fragments_sequence(
            precursors=subdict["MS1"],
            fragment_dict=subdict["MS2"],
            peak_list=subdict["peak_list"],
            ripper_dict=ripper_dict,
            err=extractor_parameters["error"],
            err_abs=extractor_parameters["err_abs"],
            ms_level=2,
            min_annotated_peak_assignment=extractor_parameters[
                "min_ms2_peak_abundance"]
        )
        
        # if fragments have been confirmed for sequence, add to confirmed
        # fragment dictionary
        if confirmed_fragments:
            confirmed_fragment_dict[sequence] = confirmed_fragments

    # return dictionary of sequences with confirmed fragments
    return confirmed_fragment_dict

def extract_MS1(
    silico_dict,
    extractor_parameters,
    ripper_dict
):
    """
    Takes an in silico sequence dict, dictionary of extractor parameters
    (passed on from input parameters .json file), mass spec data in mzml ripper
    format, and outputs a dictionary of sequences and their corresponding
    MS1 EICs (if there are MS1 hits)

    Args:
        silico_dict (dict): dictionary of sequences and corresponding masses
        extractor_parameters (dict): dictionary of extractor parameters in
            standard input parameters format
        ripper_dict (dict): dictionary of mass spec data in mzml ripper format

    Returns:
        EIC_dict (dict): dictionary of sequences and corresponding MS1 EICs
    """

    if type(extractor_parameters) == str:
        extractor_parameters = read_parameters(
            extractor_parameters)['extractor_parameters']

    if 'extractor_parameters' in extractor_parameters:
        extractor_parameters = extractor_parameters['extractor_parameters']

    # initialise dictionary to store MS1 EICs for sequences that are present
    EIC_dict = {}

    # iterate through sequences
    for sequence in silico_dict:

        # check whether silico_dict is MS or MSMS sequence dict; as format is
        # slightly different
        if type(silico_dict[sequence]) == dict:
            masses = silico_dict[sequence]["MS1"]
        else:
            masses = silico_dict[sequence]

        # make sequence masses floats, just in case they are not already
        silico_dict[sequence] = [float(mass) for mass in masses]

        # generate an MS1 EIC for sequence and associated masses
        sequence_EIC = generate_EIC(
            ions=masses,
            ms_level=1,
            ripper_dict=ripper_dict,
            err=extractor_parameters["error"],
            err_abs=extractor_parameters["err_abs"],
            min_max_intensity=extractor_parameters["min_MS1_max_intensity"],
            min_total_intensity=extractor_parameters["min_MS1_total_intensity"]
        )

        # add sequence and EIC to EIC_dict if sequence is present
        if sequence_EIC:
            EIC_dict[sequence] = sequence_EIC
            print(f'EIC generated for {sequence}')
        else:
            print(f'{sequence} not present at sufficent abundance for EIC')

    # return dict of found sequences and their corresponding MS1 EICs
    return EIC_dict

def extract_MS2(
    silico_dict,
    extractor_parameters,
    ripper_dict,
    filter_most_intense=True,
    extract_uniques=True
):
    """
    This function takes a silico MSMS sequence dict and generates MS2 EICs
    for each sequence in the dictionary, returning a dictionary of sequences
    and corresponidng MS2 EICs

    Args:
        silico_dict (dict): full MSMS silico dict in standard format
        extractor_parameters (dict): extractor parameters dict passed on
            from input parameters .json file
        ripper_dict (dict): mass spec data dict in mzml ripper format
        filter_most_intense (bool, optional): specifies whether to only
            return the most intense MS2 EIC for target sequences.
            Defaults to True.

    Returns:
        EIC_dict: dictionary of sequences and their corresponding MS2 EICs
    """
    # initiate dict to store sequences and their MS2 EICs
    EIC_dict = {}

    # retrieve MS2 min max and min total intensities
    min_max_intensity = extractor_parameters["min_MS2_max_intensity"]
    min_total_intensity = extractor_parameters["min_MS2_total_intensity"]

    # work out whether min max and min total MS2 intensities are supplied
    # as absolute intensity thresholds or decimal fractions of MS1 values
    if min_max_intensity and min_max_intensity <= 1:
        min_max_intensity = min_max_intensity*extractor_parameters[
            "min_MS1_max_intensity"]

    if min_total_intensity and min_total_intensity <= 1:
        min_total_intensity = min_total_intensity*extractor_parameters[
            "min_MS1_total_intensity"]

    # iterate through sequences, generating MS2 EICs for each of their
    # fragments
    for sequence, sequence_info in silico_dict.items():

        # retrieve sequence fragment dict
        if "MS2" in silico_dict[sequence]:
            fragments = sequence_info["MS2"]
        else:
            fragments = sequence_info

        # initate empty list for sequence MS2 EIC
        ms2_sequence_EIC = []

        # iterate through fragment ion mass lists, generating MS2 EICs for
        # each
        for fragment, masses in fragments.items():

            # generate EIC for fragment
            if fragment in sequence_info["unique_fragments"]:
                fragment_EIC = generate_EIC(
                    ions=masses,
                    ms_level=2,
                    ripper_dict=ripper_dict,
                    err=extractor_parameters["error"],
                    err_abs=extractor_parameters["err_abs"],
                    min_max_intensity=min_max_intensity,
                    min_total_intensity=min_total_intensity,
                    precursors=sequence_info["MS1"],
                    mapi=extractor_parameters["min_ms2_peak_abundance"],
                    peak_list=sequence_info["peak_list"]
                )

                if not fragment_EIC:
                    raise Exception(f'no EIC found for unique {fragment} for {sequence}')

            else:
                fragment_EIC = []

            # calculate intensity of fragment MS2 EIC
            fragment_intensity = sum(
                [Rt_I[1] for Rt_I in fragment_EIC]
            )

            # check whether only most intense EIC is to be returned
            # if not, append fragment_EIC to ms2_sequence_EIC
            if not filter_most_intense:
                ms2_sequence_EIC.append(fragment_EIC)

            # if only most intense EIC is to be returned, check whether
            # fragment_EIC is mote intense than ms2_sequence_EIC and
            # if so, reset ms2_sequence_EIC to fragment_EIC
            else:
                if fragment_intensity > sum(
                    Rt_I[1] for Rt_I in ms2_sequence_EIC
                ):
                    ms2_sequence_EIC = fragment_EIC

        # add sequence and its ms2_sequence_EIC to EIC_dict
        if ms2_sequence_EIC:
            EIC_dict[sequence] = ms2_sequence_EIC

    # return dict of sequences and ms2 EICs
    return EIC_dict
