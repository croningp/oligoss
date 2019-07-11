"""Module containing helper functions for the sequence_analyser

Reading/Writing JSON and functions for processing MSn data

.. moduleauthor:: Graham Keenan 2019
.. moduleauthor:: David Doran 2019

"""

import json

# Error thresholds
PPM = 1*(10**-6)
DEFAULT_ERROR = 10*PPM
BRUKER_ERROR = 0.01


def read_json(filename: str) -> dict:
    """Reads a JSON file and returns the data

    Arguments:
        filename {str} -- Path to the JSON file

    Returns:
        dict -- JSON data
    """

    with open(filename) as f_d:
        return json.load(f_d)


def write_json(data: dict, filename: str):
    """Writes data to a JSON file

    Arguments:
        data {dict} -- Data to write to JSON
        filename {str} -- Name of the file
    """

    with open(filename, "w") as f_d:
        json.dump(data, f_d)


def find_target_matches(
        target: float,
        mass_list: list,
        bruker: bool = True
    ) -> (list, int):
    """Performs a filter search of a mass list for a specific target

    Selects all masses in the list that are within the error thresholding

    Arguments:
        target {float} -- Mass to look for in the mass list
        mass_list {list} -- List of masses to search through

    Keyword Arguments:
        bruker {bool} -- Use Bruker thresholding or not (default: {True})

    Returns:
        list, int -- List of confirmed matches, 0 if none
    """

    # Use Bruker error if set, else default error threshold
    error = BRUKER_ERROR if bruker else DEFAULT_ERROR

    # Filter for matches that are within the threshold
    matches = list(
        filter(
            lambda x: x >= (target-error) and x <= (target+error),
            mass_list
        )
    )

    # Return 0 on no matches
    if not matches:
        return 0

    return matches


def make_EIC(matches: list, spectrum: dict) -> tuple:
    """Creates an EIC of the match intensities and retention time for that spectrum

    Totals all intensities of the matches in the spectrum, or uses 0 if no matches

    Arguments:
        matches {list} -- List of masses in the spectrum
        spectrum {dict} -- Spectrum information

    Returns:
        tuple -- EIC (retention_time, total_intensities)
    """

    ret_time = spectrum["retention_time"]

    # EIC of 0 intensity on no matches
    if not matches:
        return (ret_time, 0)

    # Total intensities of all matches
    intensity = 0
    for mass in matches:
        intensity += spectrum[f"{mass:.4f}"]

    return (ret_time, intensity)


def maximum_intensity_in_eics(eic_list: list) -> int:
    """Finds the maximum intensity in a list of lists of EICs

    Arguments:
        eic_list {list} -- List of lists containing EICs

    Returns:
        int -- Maximum intensity
    """

    # Create list of max intensities for each sublist
    maxima = [
        max([elem[1] for elem in eic]) for eic in eic_list
    ]

    # Return the maximum of that list
    return max(maxima)


def combine_duplicate_retention_times(eic_list: list) -> list:
    """Combines duplicate retention times in an EIC list and merges their intensities

    Converts the list to a dict and totals the intensities if duplicate keys exist

    Arguments:
        eic_list {list} -- EIC list

    Returns:
        list -- EIC list without duplicate intensities
    """

    out = {}

    # Go through each EIC in list
    for eic in eic_list:
        if not eic[0] in out.keys():
            # New Retention time, set value to intensity
            out[eic[0]] = eic[1]
        else:
            # Duplicate intensity, add intensity to current value
            out[eic[0]] += eic[1]

    # Convert back to list
    return [(key, value) for (key, value) in out.items()]


def check_ms2_parent_in_ms1(
        parent: str,
        ms1_masses: list,
        bruker: bool = True
    ) -> bool:
    """Checks if a MS2 parent is present in MS1 masses

    Arguments:
        parent {str} -- Parent value
        ms1_masses {list} -- List of MS1 masses

    Returns:
        bool -- Present or not
    """

    # Convert to float and search in mass list
    parent = float(parent)
    matches = find_target_matches(parent, ms1_masses, bruker=bruker)

    # Not present
    if not matches:
        return False

    return True


def get_ms2_spectrum_max_intensity(spectrum: dict) -> int:
    """Gets the maximum intensity in an MS2 spectrum

    Arguments:
        spectrum {dict} -- MS2 spectrum data

    Returns:
        int -- Maximum intensity
    """

    # Get the maximum intensity of this spectrum
    mass_list = spectrum["mass_list"]
    most_intense = max([spectrum[f"{mass:.4f}"] for mass in mass_list])

    return most_intense


def get_most_intense_ms2_peak_match(
        peak_list: list,
        spectrum: dict,
        bruker: bool = True
    ) -> int:
    """Gets the most intense intensity in a list of peak matches

    Arguments:
        peak_list {list} -- List of peaks to look for
        spectrum {dict} -- MS2 spectrum information

    Returns:
        int -- Maximum intensity
    """

    mass_list = spectrum["mass_list"]

    # Go through each peak and get maximum intensity
    most_intense = 0
    for peak in peak_list:
        matches = find_target_matches(peak, mass_list, bruker=bruker)
        if not matches:
            continue
        max_intensity = max([spectrum[f"{mass:.4f}"] for mass in matches])
        if max_intensity > most_intense:
            most_intense = max_intensity

    return most_intense


def check_ms2_mapi(
        peak_list: list,
        spectrum: dict,
        mapi_threshold: int = 90,
        bruker: bool = True
    ) -> bool:
    """Checks the Maximum Annotated Peak Intensity (MAPI) of the given spectrum

    MAPI consists of the maximum intesnity of the spectrum and maximum intensity of the peak list
    If the MAPI exceeds the threshold, succeed
    Else, fail

    Arguments:
        peak_list {list} -- List of peaks to search for
        spectrum {dict} -- MS2 spectrum information

    Keyword Arguments:
        mapi_threshold {int} -- MAPI threshold percentage (default: {90})

    Returns:
        bool -- Exceeds threshold or not
    """

    # Don't process if we have no masses
    if not spectrum["mass_list"]:
        return False

    # Get max intensity in spectrum and in peak list
    most_intense_in_spectrum = get_ms2_spectrum_max_intensity(spectrum)
    most_intense_peak = get_most_intense_ms2_peak_match(
        peak_list, spectrum, bruker=bruker
    )

    # No peak matches, just return
    if not most_intense_peak:
        return False

    # Calculate percentage of MAPI
    mapi = (most_intense_peak/most_intense_in_spectrum) * 100

    # Succeed if above threshold
    if mapi >= mapi_threshold:
        return True

    return False


def check_fragment_in_ms2(
        fragment: str,
        fragment_masses: list,
        ms2_masses: list,
        bruker: bool = True
    ) -> str:
    """Searches for a fragment in a list os MS2 masses

    Arguments:
        fragment {str} -- Fragment label
        fragment_masses {list} -- Masses for the given fragment
        ms2_masses {list} -- MS2 masses

    Returns:
        str -- fragment label if found, none otherwise
    """

    # Go through each mass for the fragment
    for mass in fragment_masses:
        # Find matches and return fragment if any found
        matches = find_target_matches(mass, ms2_masses, bruker=bruker)
        if matches:
            return fragment

    return ""


def get_combined_unique_EIC(
        unique_fragments: list,
        fragment_info: dict,
        spectra: list,
        bruker: bool = True
    ) -> list:
    """Gets a combined EIC for the unique fragments

    Searches through the spectra for unique fragments and creates a combined EIC if any found
    The most intense combined EIC is selected from the list and returned.

    Arguments:
        unique_fragments {list} -- List of unique fragment labels
        fragment_info {dict} -- Fragment information
        spectra {list} -- List of spectra

    Returns:
        list -- Most intense combined EIC of the unique fragments
    """

    combined_EIC = []

    # Go through each unique fragment
    for fragment in unique_fragments:
        EIC = []

        # Go through each spectrum
        for spectrum in spectra:
            ms2_masses = spectrum["mass_list"]
            ret_time = spectrum["retention_time"]
            fragment_masses = fragment_info[fragment]

            # Total intensities of each fragment mass if found
            intensity = 0
            for mass in fragment_masses:
                matches = find_target_matches(mass, ms2_masses, bruker=bruker)
                if not matches:
                    continue
                intensity += sum([spectrum[f"{m:.4f}"] for m in matches])

            # Construct EIC
            EIC.append((ret_time, intensity))

        # Remove duplicate retention times and make combined EIC
        EIC = combine_duplicate_retention_times(EIC)
        combined_EIC.append(EIC)

    # Select EIC based on most intense
    for pos, eics in enumerate(combined_EIC):
        combined_EIC[pos] = [eic for eic in eics if eic[1] > 0]

    return select_unique_eic(combined_EIC)


def select_unique_eic(eic_list: list) -> list:
    """Selects the most intense combined EIC

    If there's only one, return that

    Arguments:
        eic_list {list} -- List of combined EICs

    Returns:
        list -- Most intense EIC
    """

    # Return EIC if only one present
    if len(eic_list) < 2:
        return eic_list[0]

    # Get the most intense EICs and position in super list
    maxima = [
        max([(pos, elem[1]) for elem in eic])
        for (pos, eic) in enumerate(eic_list)
    ]

    # Get the position of the most intense EIC
    most_intense = max(maxima, key=lambda x: x[1])[0]

    return eic_list[most_intense]
