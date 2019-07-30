"""Module for filtering spectra based on a set of criteria
Criteria are:
    1.) Retention Time
    2.) Parent Ion (MS2 only)
    3.) Signature ions
    4.) Total intensity
    5.) Maximum intensity

.. moduleauthor:: Graham Keenan 2019
"""
def find_target(
    target,
    candidates,
    err,
    err_abs=True
):
    """
    Takes a target, associated error (absoulte mass units or ppm), list of
    candidates and returns candidates that match target within specified error

    Args:
        target (float): target m/z
        candidates (list): list of candidates for matching to target
        err (float): error threshold, in either absolute mass units or ppm
        err_abs (bool, optional): specifies whether error threshold is absolute
                        or in ppm; if True, error is in absolute mass units.
                        Defaults to True.

    Returns:
        matches (list): list of candidates that match target within specified
                        error
    """
    if not err_abs:
        err = (target*err)*10**-6

    min_hit, max_hit = target-err, target+err

    matches = filter(
        lambda x: x >= min_hit and x <= max_hit,
        candidates
    )
    return matches

def filter_retention_time(
        msdata: dict,
        ret_time_range: list
    ) -> dict:
    """Filters a spectra collection based on retention time

    Args:
        msdata (dict): Spectra collection
        ret_time_range (List[float]): Retention time range to filter for

    Returns:
        dict: Spectra that pass the filter
    """

    # List for filtered spectra
    filtered = {}

    # Iterate through each spectrum
    for key, spectrum in msdata.items():
        # Check if the retention time is within the range
        ret_time = float(spectrum["ret_time"])
        if ret_time >= ret_time_range[0] and ret_time <= ret_time_range[-1]:
            filtered[key] = spectrum

    # Return spectra sorted by retention time
    return filtered


def filter_parent_ion(
        msdata: dict,
        parent_ions: list,
        err,
        err_abs=True
    ) -> dict:
    """Filters a spectra collection based on their parent ion

    Args:
        msdata (dict): Spectra collection
        parent_ions (List[float]): List of parent ions to filter for

    Returns:
        dict: Spectra which pass the filter
    """

    # dict for filtered spectra
    filtered = {}

    # Iterate through each spectrum
    for key, spectrum in msdata.items():
        # Get the parent and all target mass matches
        parent = float(spectrum["parent"])
        ions = [float(ion) for ion in spectrum["mass_list"]]
        found_parent = find_target(parent, ions, err, err_abs)

        # If any matches are found, add the spectrum
        if found_parent:
            filtered[key] = spectrum

    # Return spectra sorted by parent
    return filtered


def filter_signature_ions(
        msdata: dict,
        sig_ions: list,
        err,
        err_abs=True
    ) -> dict:
    """Filters a spectra collection based on the presence of signature ions
    Only one signature may be found to pass the filter

    Args:
        msdata (dict): Spectra collection
        sig_ions (List[float]): List of ions to filter for
        target_func (Callable): Target find function

    Returns:
        dict: Spectra which pass the filter
    """

    # List for filtered spectra
    filtered = {}

    # Iterate through each spectrum
    for key, spectrum in msdata.items():
        # Get the mass list
        mass_list = [float(mass) for mass in spectrum["mass_list"]]

        # Check if signature ions are present
        found_ions = [
            ion for ion in sig_ions if find_target(
                ion,
                mass_list,
                err,
                err_abs)
        ]

        # If any matches are found, add to list
        if found_ions:
            filtered[key] = spectrum

    # Return spectra
    return filtered


def filter_total_intensity(
        msdata: dict,
        intensity_threshold: int
    ) -> dict:
    """Filters a spectra collection based on a total intensity threshold
    The total intensity of the spectra must exceed the threshold

    Args:
        msdata (dict): Spectra collection
        intensity_threshold (int): Total intensity threshold

    Returns:
        dict: Spectra which pass the filter
    """

    # List for filtered spectra
    filtered = {}

    # Iterate through each spectrum
    for key, spectrum in msdata.items():
        # Convert mass list to 4pt. float strings
        mass_list = [f"{m:.4f}" for m in spectrum["mass_list"]]

        # Get the total intensity of the enbtire spectrum
        total = sum([spectrum[m] for m in mass_list])

        # Exceeds the threshold, add to list
        if total >= intensity_threshold:
            filtered[key] = spectrum

    # Return spectra
    return filtered


def filter_max_intensity(
        msdata: dict,
        intensity_threshold: int
    ) -> dict:
    """Filters a spectra collection based on a maximum intensity threshold
    A spectrum's maximum intensity must exceed the threshold

    Args:
        msdata (dict): Spectra collection
        intensity_threshold (int): Maximum intensity threshold

    Returns:
        dict: Spectra which pass the filter
    """

    # List for filtered spectra
    filtered = {}

    # Iterate through each spectrum
    for key, spectrum in msdata.items():
        # Convert mass list to 4pt float strings
        mass_list = [f"{m:.4f}" for m in spectrum["mass_list"]]

        # Get the maximum intensity of the spectrum
        max_intensity = max([spectrum[m] for m in mass_list])

        # Exceeds the threshold, add to list
        if max_intensity >= intensity_threshold:
            filtered[key] = spectrum

    # Return spectra
    return filtered


def filter_mass_difference(
        msdata: dict,
        monomer_massdiffs: dict,
        total_comparisons: int,
        err,
        err_abs=True
    ) -> dict:
    """Filters the spectra based on the mass difference between peaks
    Scans through a list of monomer mass differences and check for peaks that
    match the difference.

    Args:
        msdata (dict): Spectrum data
        monomer_massdiffs (dict): Monomers and the associated mass differences
        total_comparisons (int): Total number of comparisons to perform
        target_func (Callable): Find target function

    Returns:
        dict: Spectra IDs with found monomers
    """

    # Create output map
    spectra = {}

    # Iterate through each spectrum
    for spec_id, spectrum in msdata.items():
        # Map empty list to output spectra id
        spectra[spec_id] = []

        # Get the mass list from the spectrum
        mass_list = spectrum["mass_list"]

        # Get spectrum peaks, sorted form the most intense
        spectrum_peaks = _sort_spectrum_peaks_by_intensity(spectrum)

        # Go through each monomer in the mass diffs
        for monomer, massdiff in monomer_massdiffs.items():
            # Check every peak in the spectrum, up to total comparisons
            for peak in spectrum_peaks[:total_comparisons]:
                # Go through each mass difference product
                for mdiff in massdiff:
                    # Get the peak plus/minus the mass difference
                    minus, plus = peak - mdiff, peak + mdiff

                    # Check if minus product found
                    match1 = find_target(
                        minus,
                        mass_list,
                        err,
                        err_abs)

                    # Found, log monomer and move onto next peak
                    if match1:
                        spectra[spec_id].append(monomer)
                        break

                    # No minus product, check plus product
                    match2 = find_target(
                        plus,
                        mass_list,
                        err,
                        err_abs)

                    # Found, log monomer and move onto next peak
                    if match1 or match2:
                        spectra[spec_id].append(monomer)

    # Return results
    return spectra


def _sort_spectrum_peaks_by_intensity(spectrum: dict) -> list:
    """Sorts a spectrum's peaks by the most intense peak

    Args:
        spectrum (dict): Spectra information

    Returns:
        List[float]: List of peaks sorted by the most intense
    """

    # Get the mass list from the spectrum and convert to 4pt float
    mass_list = [f"{m:.4f}" for m in spectrum["mass_list"]]

    # Map the peaks to their appropriate intensities
    peaks = {m: spectrum[m] for m in mass_list}

    # Return the peak list, sorted by the most intense
    return [
        float(peak) for (peak, _) in sorted(
            peaks.items(), key=lambda x: x[1], reverse=True
        )
    ]
