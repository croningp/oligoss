"""
Module for filtering spectra based on the following criteria:

(1) Retention Time
(2) Parent Ion (MS2 only)
(3) Signature ions
(4) Total intensity
(5) Maximum intensity

.. moduleauthor:: Graham Keenan 2019
.. moduleauthor:: David Doran 2019
"""
import time

def find_target(
    target,
    candidates,
    err,
    err_abs=True
):
    """
    Takes a target, associated error (absoulte mass units or ppm), list of
    candidates and returns candidates that match target within specified error.

    Args:
        target (float): target m/z.
        candidates (list): list of candidates for matching to target.
        err (float): error threshold, in either absolute mass units or ppm.
        err_abs (bool, optional): specifies whether error threshold is absolute
                        or in ppm; if True, error is in absolute mass units.
                        Defaults to True.

    Returns:
        matches (list): list of candidates that match target within specified
                        error.
    """

    # ensure everything is a float
    target = float(target)
    candidates = [float(candidate) for candidate in candidates]
    err = float(err)

    # check whether err is in absolute mass units or ppm
    if not err_abs:

        # convert err to ppm (if not absolute mass units)
        err = (target*err)*10**-6
        
    # calculate minimum and maximum range for matching target to candidate
    min_hit, max_hit = target-err, target+err

    # find candidates that match target within error threshold
    matches = list(filter(
        lambda x: x >= min_hit and x <= max_hit,
        candidates))

    # finally, return candidates that match target within error threshold
    return matches

def find_multiple_targets(
    targets,
    candidates,
    err,
    err_abs
):
    """
    Takes multiple target masses (m/z values) and returns matches from a list
    of candidate masses.

    Args:
        targets (list of floats): list of target m/z values.
        candidates (list of floats): list of candidate m/z values to match to
                    targets.
        err (float): error tolerance for matching targets to target masses;
                    units can either be in absolute mass units or ppm.
        err_abs (bool): specifies whether err units are absolute mass units
                    or ppm; if True, units = absolute mass units.

    Returns:
        matches (list of floats): list of candidate masses that match one or
                    more targets within error threshold specified.
    """

    # initiate list to store candidate masses that match target(s)
    matches = []
    
    # convert targets to list if not already
    if type(targets) != list:
        targets = [targets]

    # iterate through targets and return candidate masses that match target
    # masses within specified error threshold
    for target in targets:
        matches.extend(find_target(
            target,
            candidates,
            err,
            err_abs
        ))

    # return list of candidates that match ONE OR MORE target within error
    # threshold specified

    return matches

def sum_intensity_targets_spectrum(
    spectrum_dict,
    found_targets
):
    """
    This function takes a spectrum, a list of ions confirmed to be present
    in that spectrum, and returns the combined intensity / abundance of those
    ions in the spectrum.

    Args:
        spectrum_dict (dict): single spectrum dict in mzml ripper format.
        found_targets (list of floats): list of ion m/z values that have been
            found in the spectrum.

    Returns:
        intensity (float): combined intensity / abundance of all found_targets.
    """
    intensity = 0

    # ion m/z keys in spectrum_dict are strings of four point floats
    # floats may have been trimmed to one or two decimal places; so will not
    # be found in dict - this is accounted for and fixed if necessary
    for found_target in found_targets:
        try:
            intensity += float(spectrum_dict[f'{found_target}'])
        except KeyError:
            found_target = str(found_target) + '0000'
            found_target = found_target[0:found_target.find('.')+5]
            intensity += float(spectrum_dict[f'{found_target}'])

    return intensity

def filter_retention_time(
        msdata: dict,
        ret_time_range: list
    ) -> dict:
    """ Filters a spectra collection based on retention time.

    Args:
        msdata (dict): Spectra collection.
        ret_time_range (List[float]): Retention time range to filter for.

    Returns:
        dict: Spectra that pass the filter.
    """

    # iterate through each ms level and corresponding mass spec data in
    # ripper_dict
    for ms, ms_subdict in msdata.items():

        # dict for filtered spectra
        filtered = {}

        # Iterate through each spectrum
        for key, spectrum in ms_subdict.items():

            # Check if the retention time is within the range
            ret_time = float(spectrum["retention_time"])
            if ret_time >= ret_time_range[0] and ret_time <= ret_time_range[-1]:
                filtered[key] = spectrum

        msdata[ms] = filtered

    # Return spectra sorted by retention time
    return msdata

def filter_parent_ion(
        msdata: dict,
        parent_ions: list,
        err,
        err_abs=True,
        ms_level=2
    ) -> dict:
    """ Filters a spectra collection based on their parent ion.

    Args:
        msdata (dict): Spectra collection.
        parent_ions (List[float]): List of parent ions to filter for.

    Keyword Args:
        ms_level (int): specifies ms_level for spectra being screened (default:
                    2).
    Returns:
        dict: Spectra which pass the filter.
    """

    if ms_level < 2:
        raise Exception(f'ms_level for parent ion filter must be 2 or greater, current ms_level set to {ms_level}')

    # initialise dict to store filtered spectra, keeping all spectra at
    # ms levels not selected for filtering
    filtered_dict = {
        ms : msdata[ms]
        for ms in msdata
        if ms != ms_level
    }
    filtered_dict.update({f'ms{ms_level}': {} })

    # Iterate through each spectrum
    for key, spectrum in msdata[f'ms{ms_level}'].items():

        # Get the parent and all target mass matches
        parent = float(spectrum["parent"])

        found_parents = find_multiple_targets(
            targets=parent_ions,
            candidates=[parent],
            err=err,
            err_abs=err_abs
        )

        # If any matches are found, add the spectrum
        if found_parents:
            filtered_dict[f'ms{ms_level}'][key] = spectrum

    # if a so much as smell the hints of a deep copy, al slit yer throat
    # Return spectra sorted by parent
    return filtered_dict

def filter_signature_ions(
        msdata: dict,
        sig_ions: list,
        err,
        err_abs=True,
        ms_level=2
    ) -> dict:
    """ Filters a spectra collection based on the presence of signature ions
    Only one signature may be found to pass the filter.

    Args:
        msdata (dict): Spectra collection.
        sig_ions (List[float]): List of ions to filter for.
        target_func (Callable): Target find function.

    Keyword Args:
        ms_level (int): specifies MS level for signature ions; (default = 2).
    Returns:
        dict: Spectra which pass the filter.
    """

    # List for filtered spectra
    filtered = {}

    msdict = msdata[f'ms{ms_level}']

    # Iterate through each spectrum
    for key, spectrum in msdict.items():
        # Get the mass list
        mass_list = [float(mass) for mass in spectrum["mass_list"]]

        # Check if signature ions are present
        found_ions = [
            ion for ion in sig_ions if find_target(
                target=ion,
                candidates=mass_list,
                err=err,
                err_abs=err_abs)
        ]

        # If any matches are found, add to list
        if found_ions:
            filtered[key] = spectrum

    msdata[f'ms{ms_level}'] = filtered

    # Return msdata with spectra not containing signature ions at specified
    # ms level removed
    return msdata

def filter_total_intensity(
        msdata: dict,
        intensity_threshold: float,
        ms_level: int
    ) -> dict:
    """ Filters a spectra collection based on a total intensity threshold
    The total intensity of the spectra must exceed the threshold.

    Args:
        msdata (dict): Spectra collection.
        intensity_threshold (int): Total intensity threshold.

    Returns:
        dict: Spectra which pass the filter.
    """
    print(f'total intensity threshold for ms{ms_level} = {intensity_threshold}')
    for ms, msdict in msdata.items():

        # List for filtered spectra
        filtered = {}

        if int(ms[-1]) != ms_level:
            filtered = msdict

        elif int(ms[-1]) == ms_level:
            # Iterate through each spectrum
            for key, spectrum in msdict.items():
                # Convert mass list to 4pt. float strings
                mass_list = [f"{m:.4f}" for m in spectrum["mass_list"]]

                # Get the total intensity of the enbtire spectrum
                total = sum([spectrum[m] for m in mass_list])

                # Exceeds the threshold, add to list
                if total >= intensity_threshold:
                    filtered[key] = spectrum

        msdata[f'ms{ms_level}'] = filtered

    # Return spectra
    return msdata

def filter_max_intensity(
        msdata: dict,
        intensity_threshold: float,
        ms_level: int
    ) -> dict:
    """ Filters a spectra collection based on a maximum intensity threshold
    A spectrum's maximum intensity must exceed the threshold.

    Args:
        msdata (dict): Spectra collection.
        intensity_threshold (int): Maximum intensity threshold.

    Returns:
        dict: Spectra which pass the filter.
    """
    print(f'max intensity filter for ms{ms_level} = {intensity_threshold}')
    for ms, msdict in msdata.items():

        # dict for filtered spectra
        filtered = {}

        # check ms level of spectra
        if int(ms[-1]) != ms_level:
            filtered = msdict
        elif int(ms[-1]) == ms_level:
            # Iterate through each spectrum
            for key, spectrum in msdict.items():
                # Convert mass list to 4pt float strings
                mass_list = [f"{m:.4f}" for m in spectrum["mass_list"]]

                # Get the maximum intensity of the spectrum
                max_intensity = max([spectrum[m] for m in mass_list])

                # Exceeds the threshold, add to list
                if max_intensity >= intensity_threshold:
                    filtered[key] = spectrum

        msdata[ms] = filtered

    # Return spectra
    return msdata

def find_ms2_signature_ions(
    monomer_list,
    signature_ion_type, 
    spectrum, 
    error, 
    error_abs, 
    signature_ion_dict
    ):
    """ This function searches an MS2 spectrum and looks for signature ions
    from the monomer list. If a dominant signature ion is found, the monomer is 
    confirmed.

    Arguments:
        monomer_list {list} -- list of monomers within the sample.

    Returns:
        list - list of monomers with confirmed signature ions
    """
    # initiate confirmed monomer dict
    confirmed_monomers = []

    # keep spectra only (remove retention time, mass lists, parent peak)
    keys_to_remove = ["retention_time", "mass_list", "parent"]
    spectra_only = {k:v for k,v in spectrum.items() if k not in keys_to_remove}

    for monomer in monomer_list:

        # get mass for dominant signature ion
        for signature_ions in signature_ion_dict[signature_ion_type]:
            if signature_ions[0] == monomer:
                signature_frag_mass = signature_ions[1]

                # check if signature fragment mass is in the spectrum
                # if so, add to dictionary
                signature_frag_search = find_multiple_targets(
                    signature_frag_mass,
                    candidates=spectra_only.keys(),
                    err=error,
                    err_abs=error_abs
                )

                if signature_frag_search:
                    confirmed_monomers.append(monomer)

    # return dictionary of confirmed monomers, their confirmed fragment type 
    # and list of unconfirmed monomers
    return confirmed_monomers


def apply_dynamic_exclusion_peaks(
    peak_list=[[200,1], [202, 2], [201, 1], [100,1], [101, 2]],
    exclusion_window=1,
    exclusion_units_abs=False
):

    # ensure peak list is sorted by intensity, in descending order - with most
    # intense peak first in the list 
    peak_list = sorted(
        peak_list, 
        key=lambda peak: peak[1],
        reverse=True
    )
    
    most_intense = peak_list[0]

    filtered_peaks = most_intense

    range_match = []

    while not range_match:
        
        excluded_peaks = find_multiple_targets(
            targets=[most_intense[0]],
            candidates=[peak[0] for peak in peak_list[1::]],
            err=exclusion_window,
            err_abs=exclusion_units_abs)

        if excluded_peaks: 
            raise Exception(f'excluded peaks = {excluded_peaks}')
        trimmed_peak_list = [
            peak for peak in peak_list if peak[0] not in excluded_peaks]
        
        filtered_peaks.extend(trimmed_peak_list)
        
        peak_list = trimmed_peak_list

        most_intense, least_intense = peak_list[0], peak_list[-1]

        print(f'peak_list = {peak_list}')

        range_match = find_target(
            target=most_intense[0],
            candidates=[least_intense[0]],
            err=exclusion_window,
            err_abs=exclusion_units_abs)
        
        
    filtered_peaks = sorted(list(
        set(filtered_peaks), key = lambda peak: peak[1], reverse=True))
    
    return filtered_peaks
        
def find_common_peaks_massdiff_spectra(spectral_assignments: dict) -> dict:
    """
    This function takes a dict of spectral assignments, with spectra identified
    as containing monomer-specific mass differences and / or signature ions and
    finds peaks that match potential ladders for two or more monomers. 
    NOTE: for full format of this dict, see output of 
    'filter_monomer_fingerprints' function (also in extractor_helpers.py). 
    
    Args:
        spectral_assignments (dict): dict of spectrum ids and associated 
            monomer massdiffs and / or signatures that have been identified
            from 'filter_monomer_fingerprints' function. Format: 
            {
                'spectrum_id': {
                    'retention_time': retention_time (float),
                    'precursor': precursor (float),
                    'monomer1': {
                        'signatures': [peak (float),...]
                        'mass_diffs': [peak (float),...]
                    },
                    'monomer2': {
                        'signatures': [peak (float),...],
                        'mass_diffs': [peak (float),...]
                    },
                    ... 
                },
                ...
            }
            where signatures list contains specific signature ions identified 
            for monomer and mass_diffs list contains peaks which are separated
            from one or more peak by a mass corresponding to target monomer. 
    
    Returns:
        dict: dict of spectral assignments, with an extra key-value pair added
            for peaks in 'mass_diffs' list which are common to two or more 
            monomers
    """

    # iterate through spectra and their monomer fingerprint assignments 
    for spectrum_id, mdiffs in spectral_assignments.items():
        
        # identify mass_diff peaks for eahc monomer 
        non_monomer_keys = [
            'retention_time', 'precursor', 'base_peak', 'base_peak_assignment']
        print(f'mdiffs = {mdiffs}')

        monomer_massdiffs = {
            monomer: sorted(list(set(mdiffs[monomer]['mass_diffs'])), reverse=True)
            for monomer in mdiffs
            if monomer not in non_monomer_keys}

        # init lists to store ALL identified peaks for all monomers and peaks 
        # that are common to two or more monomers
        all_peaks, common_peaks = [], []
        
        # check whether info for more than one monomer has been found; if not,
        # no common peaks will be found so pass
        if len([key for key in monomer_massdiffs]) == 1:
            pass 
        else:

            # iterate through mass diff peaks for each monomer and add to list
            # of ALL peaks
            for peaks in monomer_massdiffs.values():
                all_peaks.extend(peaks)
            
            # for each peak, check whether it occurs more than once; if so, it
            # must be found in mass diffs for more than one monomer - add to
            # common peaks list 
            for peak in all_peaks: 
                if all_peaks.count(peak) > 1: 
                    common_peaks.append(peak)
        
        # remove duplicates from common peaks and sort in descending order 
        common_peaks = sorted(list(set(common_peaks)), reverse=True)

        # update input dict with 'common_peaks' and list of common peaks
        spectral_assignments[spectrum_id].update(
            {'common_peaks': common_peaks})

    return spectral_assignments 

def _sort_spectrum_peaks_by_intensity(spectrum: dict) -> list:
    """ Sorts a spectrum's peaks by the most intense peak.

    Args:
        spectrum (dict): Spectra information.

    Returns:
        List[list]: List of peaks sorted by the most intense. NOTE: peaks are
            in format [mz, I] where mz, I = m/z and intensity value of peak

    """

    # Get the mass list from the spectrum and convert to 4pt float
    mass_list = [f"{m:.4f}" for m in spectrum["mass_list"]]

    # Map the peaks to their appropriate intensities
    peaks = {m: spectrum[m] for m in mass_list}

    # return list of peaks, sorted (in descending order) by intensity
    return sorted(
        [peak for peak in peaks.items()],
        key = lambda peak: float(peak[1]),
        reverse=True)

def return_most_intense_peaks_spectrum(
    spectrum,
    N_peaks=None,
    min_relative_intensity=0.01,
    dynamic_exclusion_window=None,
    dynamic_exclusion_abs=False
):
    """
    Takes a spectrum and returns N most intense peaks as a list of sublists 
    in format [[mz, I]...] where mz = m/z of peak, I = intensity (abundance)
    
    Args:
        spectrum (dict): single spectrum dict in mzml ripper format
        N_peaks (int): number of most intense peaks to return
    
    Returns:
        list: list of most intense peaks, in descending order of intensity 
    """
    spectrum_peaks = _sort_spectrum_peaks_by_intensity(spectrum=spectrum)
    min_abs_intensity = spectrum_peaks[0][1]*min_relative_intensity 

    spectrum_peaks = [
        peak for peak in spectrum_peaks
        if peak[1] >= min_abs_intensity]

    if dynamic_exclusion_window:
        spectrum_peaks = apply_dynamic_exclusion_peaks(
            peak_list=spectrum_peaks,
            exclusion_window=dynamic_exclusion_window,
            exclusion_units_abs=dynamic_exclusion_abs)

    if N_peaks: 
        spectrum_peaks = spectrum_peaks[:min([N_peaks, len(spectrum_peaks)])]
    
    return spectrum_peaks
    

def apply_pre_screening_filters(
    ripper_dict,
    min_rt,
    max_rt,
    essential_signatures,
    signature_ms_level,
    precursor_ions,
    precursor_ion_ms_level,
    min_MS1_total_intensity,
    min_MS2_total_intensity,
    min_MS1_max_intensity,
    min_MS2_max_intensity,
    err,
    err_abs
):
    """
    Takes an MS dictionary that is read from an mzml ripper json and
    returns a filtered dictionary of spectra within the input ripper_dict that
    meet the filter criteria specified (described below).

    Args:
        ripper_dict (dict): mass spec data dict in mzml ripper format.
        min_rt (float): minimum retention time for spectra (in minutes).
        max_rt (float): maximum retention time for spectra (in minutes).
        essential_signatures (list): list of monomers that MUST have signature
            ions present in all spectra screened. Only spectra containing
            these monomer signature ions at the specified ms level will be
            returned.
        signature_ms_level (int): ms level for monomer signature screening.
        precursor_ions (list): list of precursor ions. At least one of these
            must be the precursor for any MSn spectra.
        precursor_ion_ms_level (int): MS level for screening precursors. This
            MUST be greater than 1, as MS1 spectra have no precursors.
        min_MS1_total_intensity (float): minimum total intensity of MS1 spectra.
            Spectra with a total intensity below this value will be filtered
            out.
        min_MS2_total_intensity (float): minimum total intensity of MS2 spectra.
            Spectra with a total intensity below this value will be filtered
            out. NOTE: This is either given as an absolute value or a decimal
            fraction of min_MS1_total_intensity; either is acceptable.
        min_MS1_max_intensity (float): minimum maximum intensity of MS1 spectra.
            Spectra with a maximum intensity below this value will be filtered
            out.
        min_MS2_max_intensity (float): minimum maximum intensity of MS2 spectra.
            Spectra with a maximum intensity below this value will be filtered
            out. NOTE: This is either given as an absolute value or a decimal
            fraction of min_MS1_max_intensity; either is acceptable.
        err (float): error threshold for screening - either in absolute mass
            units or ppm.
        err_abs (bool): specifies whether err units are absolute mass units
            or ppm.

    Returns:
        dict: A filtered dictionary of spectra within the input ripper_dict that 
        meet the filter criteria specified.
    """
    # check for min, max rt and apply filter
    if not min_rt and not max_rt:
        ripper_dict = ripper_dict
        print('ripper dict not filtered by retention time')
    else:
        if not min_rt:
            min_rt = 0
        if not max_rt:
            max_rt = math.inf
        ripper_dict = filter_retention_time(
            msdata=ripper_dict,
            ret_time_range=[min_rt, max_rt]
        )

    # check for essential signatures; if any are specified, filter out
    # spectra that do not contain these signatures at specified ms level
    if essential_signatures:
        for signature in essential_signatures:
            ripper_dict = filter_signature_ions(
                msdata=ripper_dict,
                sig_ions=[signature],
                err=err,
                err_abs=err_abs,
                ms_level=signature_ms_level)

    # check for precursor ions; if any are specified, filter out any MSn
    # spectra that do not match any of these precursors
    if precursor_ions:
        ripper_dict = filter_parent_ion(
            msdata=ripper_dict,
            parent_ions=precursor_ions,
            err=err,
            err_abs=err_abs,
            ms_level=precursor_ion_ms_level
        )

    # check for minimum MS1 total intensity; if specified, filter out any MS1
    # spectra that have a total intensity below this value
    if min_MS1_total_intensity:
        print(f'filtering out MS1 spectra with total I below {min_MS1_total_intensity}')
        ripper_dict = filter_total_intensity(
            msdata=ripper_dict,
            intensity_threshold=min_MS1_total_intensity,
            ms_level=1
        )

    # check for minimum MS2 total intensity; if specified, filter out any MS2
    # spectra that have a total intensity below this value
    # also check to see if min_MS2_total_intensity is an absolute value or a
    # decimal fraction of min_MS1_total_intensity
    if min_MS2_total_intensity:
        if min_MS2_total_intensity < 1:
            min_MS2_total_intensity = min_MS2_total_intensity * min_MS2_total_intensity
        ripper_dict = filter_total_intensity(
            msdata=ripper_dict,
            intensity_threshold=min_MS2_total_intensity,
            ms_level=2
        )


    # check for minimum MS1 max intensity; if specified, filter out any MS1
    # spectra that have a total intensity below this value
    if min_MS1_max_intensity:
        ripper_dict = filter_max_intensity(
            msdata=ripper_dict,
            intensity_threshold=min_MS1_max_intensity,
            ms_level=1
        )

    # check for minimum MS2 max intensity; if specified, filter out any MS2
    # spectra that have a maximum intensity below this value
    # also check to see if min_MS2_max_intensity is an absolute value or a
    # decimal fraction of min_MS1_max_intensity
    if min_MS2_max_intensity:
        if min_MS2_max_intensity < 1:
            min_MS2_max_intensity = min_MS2_max_intensity*min_MS1_max_intensity
        ripper_dict = filter_max_intensity(
            msdata=ripper_dict,
            intensity_threshold=min_MS2_max_intensity,
            ms_level=2
        )
    return ripper_dict

def find_massdiffs_spectrum_peaks(
    spectrum_peaks,
    massdiffs,
    err,
    err_abs
):

    massdiff_matches, masses = [], [peak[0] for peak in spectrum_peaks]
    len_masses = len(masses)
    for massdiff in massdiffs: 

        for peak in spectrum_peaks: 

            high, low = float(peak[0]) + massdiff, float(peak[0]) - massdiff 
            
            high_matches = find_target(
                target=high,
                candidates=masses,
                err=err,
                err_abs=err_abs)
            
            low_matches = find_target(
                target=low,
                candidates=masses,
                err=err,
                err_abs=err_abs)
            
            if high_matches:
                massdiff_matches.append(high)
                massdiff_matches.extend(high_matches)
            
            if low_matches:
                massdiff_matches.append(low)
                massdiff_matches.extend(low_matches)
        
        if len(masses) != len_masses:
            raise Exception(f'n masses has been reduced from {len_masses} to {len(masses)}')
    
    return sorted(list(set(massdiff_matches)), reverse=True)