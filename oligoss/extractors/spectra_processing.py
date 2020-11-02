"""
This file contains functions for processing mass spectra independently of
any in silico data for target ions.
"""
import itertools

def get_mass_defect_peaks(
    peaks,
    base_unit=None
):
    """
    Takes list of peaks and returns fractional mass defect or Kendrick mass
    defect.

    Args:
        peaks (List[Tuple[float]]): list of peaks in format: [(m/z, I)...]
            where m/z, I = mass-to-charge ration, intensity.
        base_unit (float, optional): neutral monoisotopic mass of suspected
            repeating unit. NOTE: if specified, list of
            Kendrick Mass Defects will be returned; otherwise, fractional
            mass defect will be returned. Defaults to None.

    Returns:
        List[Tuple[float]]: list of mass defects in format: [(m/z, defect)...]
    """

    #  base unit has been specified, so return Kendrick Mass Defects
    if base_unit:
        return [
            (peak[0], peak[0] * (round(base_unit) / base_unit))
            for peak in peaks
        ]

    #  no base unit has been specified, so return fractional mass defects
    return [
        (peak[0], peak[0] - round(peak[0]))
        for peak in peaks
    ]


def generate_ion_chromatogram(
    spectra: list,
    chromatogram_type: str,
    exclusion_window: float,
    min_rel_intensity: float,
    min_abs_intensity: float
) -> list:
    """
    Takes a list of spectra and returns either Total Ion Current (TIC) or
    Base Peak Chromatogram (BPC).

    Args:
        spectra (List[dict]): list of ripper spectra.
        chromatogram_type (str): either "tic" or "bpc" to return TIC or BPC.
        exclusion_window (float): exclusion window for binning peaks, in units
            of parts per million (ppm).
        min_rel_intensity (float): minimum relative intensity (relative to base
            peak) for peaks to be included in binning / intensity contributions.
            NOTE: expressed as a decimal fraction between 0 and 1.
        min_abs_intensity (float): minimum absolute intensity for peaks to be
            included in binning / intensity contributions.

    Raises:
        Exception: raised if chromatogram_type not in ["tic", "bpc"].

    Returns:
        List[list]: list of retention times and associated ion info. Formats:
            1. For BPC: [[rt, m/z, I]...] where rt, m/z and I = retention time,
                mass-to-charge ratio of most intense peak @ rt, and I =
                intensity of peak in m/z bin.
            2. For TIC: [[rt, I]...] where rt, I = retention time and I = total
                intensity of all peaks @ rt.
    """

    #  init list to store final chromatogram
    chromatogram = []

    #  iterate through spectra and append to chromatogram
    for spectrum in spectra:
        chrom_peak = [spectrum["retention_time"]]

        #  sort peaks and apply exclusion window
        spectrum_peaks = bin_spectrum_peaks(
            spectrum=spectrum,
            ppm_window=exclusion_window,
            min_rel_intensity=min_rel_intensity,
            min_abs_intensity=min_abs_intensity,
            intensity_grouping="sum"
        )

        #  depending on chromatogram type, append appropriate data to final
        #  chromatogram
        if chromatogram_type.lower() == "tic":
            if spectrum_peaks:
                chrom_peak.append(sum([
                    x[1] for x in spectrum_peaks
                ]))
            else:
                chrom_peak.append(0)
        elif chromatogram_type.lower() == "bpc":
            if spectrum_peaks:
                chrom_peak.extend(
                    spectrum_peaks[0]
                )
            else:
                chrom_peak.extend([0, 0])
        else:
            raise Exception("chromatogram_type must be either 'tic' or 'eic'")

        chromatogram.append(chrom_peak)

    return chromatogram

def bin_spectrum_peaks(
    spectrum: dict,
    ppm_window: float,
    min_abs_intensity: float,
    min_rel_intensity: float,
    intensity_grouping: str
):
    """
    Takes a spectrum and groups spectrum peaks within specified ppm window.

    Args:
        spectrum (dict): spectrum dict in ripper format.
        ppm_window (float): ppm window for matching m/z values.
        min_abs_intensity (float): minimum absolute intensity for peaks to
            be considered.
        min_rel_intensity (float): minimum relative intensity for peaks to
            be considered. NOTE: this MUST be expressed as a decimal fraction
            in range 0 -> 1.
        intensity_grouping (str): specifies whether intensity of binned peaks
            are represented as mean or sum of binned intensities. NOTE: there
            are two options for this: "mean" or "sum" for mean or summed
            intensities, respectively.

    Returns:
        dict: spectrum dict with peaks grouped within ppm window.
    """
    #  get list of spectrum peaks, sorted in DESCENDING order by intensity
    peaks = list_spectrum_peaks(spectrum)

    #  bin spectrum peaks by mass, matched within ppm_window
    peaks = bin_peaks_by_mz(
        peaks=peaks,
        ppm_window=ppm_window,
        min_rel_intensity=min_rel_intensity,
        min_abs_intensity=min_abs_intensity,
        intensity_grouping=intensity_grouping)

    #  create new spectrum dict with binned peaks
    binned_spectrum = {str(peak[0]): peak[1] for peak in peaks}
    binned_spectrum["mass_list"] = [float(peak) for peak in binned_spectrum]

    #  add precursor if applicable
    if "parent" in spectrum:
        binned_spectrum["parent"] = spectrum["parent"]

    return binned_spectrum

def bin_peaks_by_mz(
    peaks: list,
    ppm_window: float,
    min_rel_intensity: float,
    min_abs_intensity: float,
    intensity_grouping: str
) -> list:
    """
    Takes a peak list and returns list of peaks binned within ppm_window.

    Args:
        peaks (List[List[int]]): list of peaks in format:
            [[m/z, I], ...] where m/z, I = mass-to-charge ratio and absolute
            intensity of peak, respectively.
        ppm_window (float): binning window in parts per million (ppm).
        min_rel_intensity (float): minimum relative intensity for peaks to be
            included in binning. Set relative to base peak. NOTE: must be
            expressed as decimal fraction between 0 and 1.
        min_abs_intensity (float): minimum absolute intensity for peaks to be
            included in binning.
        intensity_grouping (str): specifies whether final peak intensities
            are the sum of combined peaks or an average intensity. NOTE: this
            arg must be either "sum" or "mean" for summing and averaging peak
            intensities respectively.

    Returns:
        List[List[float]]: list of peaks sorted in descending order by intensity
            in format: [[m/z, I]...]
            where m/z is the mass-to-charge ratio of most intense peak in bin
            and I is cumulative intensity of all peaks within bin (i.e. peaks
            that match m/z within ppm_window).
    """

    #  sort peaks by intensity (DESCENDING order - i.e. most intense first)
    peaks = sorted(peaks, key=lambda x: x[1], reverse=True)

    #  remove peaks which fall below absolute and relative intensity thresholds
    if min_abs_intensity:
        peaks = [
            peak for peak in peaks if peak[1] >= min_abs_intensity
        ]
    if min_rel_intensity:
        peaks = [
            peak for peak in peaks if peak[1] >= peaks[0][1] * min_rel_intensity
        ]

    #  no peaks or only one peak, so nothing to bin. Returns peaks as they are
    if not peaks or len(peaks) == 1:
        return peaks

    def bin_peaks():
        """
        Inner function for binning peaks within error threshold (ppm)
        """

        #  get target peak to set bin
        most_intense = peaks[iter_count]

        #  find masses of peaks to be binned with target peak
        binned_peaks = find_target_match(
            target_mass=most_intense[0],
            candidate_masses=[x[0] for x in peaks[iter_count::]],
            error=ppm_window,
            err_units="ppm"
        )
        if len(binned_peaks) > 1:

            #  remove peaks that have been binned with target
            excluded_peaks = [
                x for x in peaks
                if x[0] in binned_peaks and x[0] != most_intense[0]]

            #  get list of intensities for peaks which match target peak
            #  within ppm_window
            matching_peaks = [x[1] for x in peaks if x[0] in binned_peaks]

            #  check the intensity grouping method and group accordingly - final
            #  peak intensities are either sum or mean of binned peaks
            if intensity_grouping.lower() == "sum":
                peaks[iter_count][1] = sum(matching_peaks)
            elif intensity_grouping.lower() == "mean":
                peaks[iter_count][1] = sum(matching_peaks) / len(matching_peaks)
            else:
                raise Exception(
                    "intensity_grouping must be either 'sum' or 'mean'")
            for peak in excluded_peaks:
                peaks.remove(peak)

        return peaks

    #  init counter to track how many peak bins have been found
    iter_count = 0

    #  peak binners gonna bin peaks
    while iter_count < len(peaks):
        peaks = bin_peaks()
        iter_count += 1

    return sorted(peaks, key=lambda x: x[1], reverse=True)

def sort_spectrum_peaks_by_intensity(
    spectrum: dict,
    min_rel_intensity: float,
    min_abs_intensity: float
) -> list:
    """
    Takes a single ripper spectrum and returns list of peaks in spectrum sorted
    in descending order by intensity.

    Args:
        spectrum (dict): single ripper spectrum dict.
        min_rel_intensity (float): minimum relative intensity of peaks. NOTE:
            this must be expressed as a decimal fraction between 0 and 1.
        min_abs_intensity (float): minimum absolute intensity of peaks.

    Raises:
        Exception: raised if min_rel_intensity > 1.

    Returns:
        List[list]: list of peaks in format [[m/z, i]...] where m/z, i = mass to
            charge ratio and intensity of peaks. Sorted in descending order by
            intensity.
    """
    #  get list of spectrum peaks, sorted in DESCENDING order by intensity
    peaks = sorted(list_spectrum_peaks(spectrum), lambda x: x[1], reverse=True)

    #  apply absolute intensity threshold
    if min_abs_intensity:
        peaks = [peak for peak in peaks if peak[1] >= min_abs_intensity]
        if not peaks:
            return []

    #  apply relative intensity threshold
    if min_rel_intensity:
        if min_rel_intensity > 1:
            raise Exception(
                "min_rel_intensity cannot be greater than 1\n"
                f"current value is {min_rel_intensity}")
        return [
            peak for peak in peaks
            if peak[1] >= peaks[0][1] * min_rel_intensity]

    return peaks

def list_spectrum_peaks(spectrum: dict) -> list:
    """
    Takes a ripper spectrum and retrieves list of peaks in spectrum.

    Args:
        spectrum (dict): spectrum in ripper dict format.
    Returns:
        List[List[float]]: list of raw peaks in spectrum in format:
            [[m/z, I], ...] where m/z, I = mass-to-charge ratio, absolute
            intensity of peaks, respectively.
    """
    #  get list of unique masses in spectrum
    spectrum_masses = [f"{m:.4f}" for m in spectrum["mass_list"]]

    #  get list of peaks in format [[mass, intensity]...] sorted in DESCENDING
    #  order by intensity
    return [[float(m), spectrum[m]] for m in spectrum_masses]

def find_target_match(
    target_mass: float,
    candidate_masses: list,
    error: float,
    err_units: str
) -> list:

    if err_units == "ppm":
        error = (target_mass / 1E6) * error
    upper, lower = target_mass + error, target_mass - error

    return [
        m for m in candidate_masses if m >= lower and m <= upper
    ]

def combine_spectra_peaks(
    spectra,
    min_rel_intensity,
    min_abs_intensity,
    ppm_window=None,
    intensity_grouping="mean"
):
    """
    Takes a list of spectra and returns list of all peaks in spectra that
    exceed intensity thresholds.

    Args:
        spectra (List[dict]): list of ripper spectra dicts.
        min_rel_intensity (float): minimum relative intensity for peaks to be
            considered. NOTE: this is expressed as a decimal fraction with
            range 0 -> 1. Peaks with a relative intensity below
            (min_rel_intensity * base_peak) will be discarded.
        min_abs_intensity (float): minimum absolute intensity for peaks to be
            considered. Peaks with an absolute intensity below this value will
            be discarded.
        ppm_window (float, optional): specifies ppm window for combining peaks.
            If None, raw peaks will be combined as they are with only filtering
            for intensity thresholds.
        intensity_grouping (str, optional): specifies whether to combine
            peaks by summing or averaging intensity. NOTE: this must be either
            "mean" or "sum" for averaging or summing peak intensities,
            respectively. This is only used if a ppm_window is specified for
            combining peaks.

    Returns:
        List[List[float]]: list of peaks in format:
            [[m/z, I], ...] where m/z, I = mass-to-charge ratio and absolute
            intensity of peak, respectively.
    """

    #  init list to store peaks
    peaks = []

    #  iterate through spectra, collecting peaks
    for spectrum in spectra:
        peaks.extend(sort_spectrum_peaks_by_intensity(
            spectrum=spectrum,
            min_rel_intensity=min_rel_intensity,
            min_abs_intensity=min_abs_intensity)
        )

    if ppm_window:
        peaks = bin_peaks_by_mz(
            peaks=peaks,
            ppm_window=ppm_window,
            min_rel_intensity=None,
            min_abs_intensity=None,
            intensity_grouping=intensity_grouping
        )

    #  return peaks, sorted in descending order by intensity
    return sorted(peaks, lambda x: x[1], reverse=True)

def make_spectrum_from_peaks(
    peaks: list,
    parent=None
) -> dict:
    """
    Takes a list of peaks and generates a spectrum in ripper format from peaks.
    Useful for manipulating spectra using subsets of peaks or modified peaks.

    Args:
        peaks (List[List[int]]): list of peaks in format:
            [[m/z, I], ...] where m/z, I = mass-to-charge ratio and absolute
            intensity for peaks, respectively.
        parent (float, optional): exact precursor / parent m/z value, if
            applicable. Not relevant for MS1 spectra. Defaults to None:float.

    Returns:
        dict: spectrum dict in ripper format.
    """

    #  generate spectrum dict
    spectrum = {
        str(peak[0]): peak[1]
        for peak in peaks}

    #  add mass list to dict
    spectrum["mass_list"] = [
        peak[0] for peak in peaks]

    #  if there is a precursor / parent mass in spectrum, add that info to
    #  spectrum dict
    if parent:
        spectrum["parent"] = parent

    return spectrum

def generate_consensus_spectrum(
    spectra: list,
    min_rel_intensity: float,
    min_abs_intensity: float,
    min_peak_identity: float,
    ppm_window: float,
    intensity_grouping: str
) -> dict:
    """
    Takes a list of spectra in ripper format and returns a consensus spectrum in
    same format.

    Args:
        spectra (List[dict]): list of ripper spectra.
        min_rel_intensity (float): minimum relative intensity for peaks to
            be considered. NOTE: must be decimal fraction in range 0 -> 1.
        min_abs_intensity (float): minimum absolute intensity for peaks to be
            considered.
        min_peak_identity (float): minimum % of spectra for peaks to be present
            in to be included in final consensus spectrum. NOTE: this must
            be expressed as a decimal fraction in range 0 -> 1.
        ppm_window (float): error window for grouping peaks by m/z (units = ppm)
        intensity_grouping (str): specifies whether final intensities of grouped
            peaks in consensus spectrum should equal the sum or mean of
            individual peak intensities. NOTE: must be either "sum" or "mean" to
            specify mean or summed intensities, respectively.

    Returns:
        dict: consensus spectrum in ripper dict format.
    """
    #  combine peak lists for all spectra
    peaks = itertools.chain([
        list_spectrum_peaks(spectrum) for spectrum in spectra])

    #  apply peak binning
    peaks = bin_peaks_by_mz(
        peaks=peaks,
        ppm_window=ppm_window,
        min_rel_intensity=min_rel_intensity,
        min_abs_intensity=min_abs_intensity,
        intensity_grouping=intensity_grouping)

    #  make sure to record precursor m/z if inputs are MS2-n spectra. This is
    #  required for returning final spectra in correct format.
    if "parent" in spectra[0]:
        parent = spectra[0]["parent"]
    else:
        parent = None

    #  if peaks must be present in minimum % of spectra, remove those that do
    #  not meet this threshold
    if min_peak_identity:
        peak_counts = {}
        min_peak_count = min_peak_identity * len(spectra)
        for peak in peaks:
            count = 0
            for spectrum in spectra:
                raw_peaks = sort_spectrum_peaks_by_intensity(
                    spectrum=spectrum,
                    min_rel_intensity=min_rel_intensity,
                    min_abs_intensity=min_abs_intensity)
                if find_target_match(
                        target_mass=peak[0],
                        candidate_masses=[x[0] for x in raw_peaks],
                        error=ppm_window,
                        err_units="ppm"):
                    count += 1
                if count >= min_peak_count:
                    peak_counts[peak[0]] = count
                    break
        return make_spectrum_from_peaks(
            peaks=[peak for peak in peaks if peak[0] in peak_counts],
            parent=parent)

    return make_spectrum_from_peaks(peaks=peaks, parent=parent)

def calculate_ma_spectrum_raw(spectrum: dict):
    """
    Calculate MA for spectrum. NOTE: if calculating MA of consensus spectra or
    intensity filters are required, pass in appropriately processed spectrum.

    Args:
        spectrum (dict): spectrum dict in ripper format

    Returns:
        float: MA value for spectrum (i.e. number of unique peaks)
    """
    return len([mz for mz in spectrum if mz not in ["parent", "mass_list"]])

def calculate_ma_consensus_spectra(
    spectra: list,
    ppm_window: float,
    intensity_grouping: str,
    min_abs_intensity: float,
    min_rel_intensity: float,
    min_peak_identity: float
) -> float:
    """
    Takes a list of spectra in ripper format, groups into a single consensus
    spectrum and then returns the calculated MA of that consensus spectrum.

    Args:
        spectra (List[dict]): list of ripper spectra.
        ppm_window (float): error window for grouping peaks by m/z
            (units = ppm).
        intensity_grouping (str): specifies whether final intensities of grouped
            peaks in consensus spectrum should equal the sum or mean of
            individual peak intensities. NOTE: must be either "sum" or "mean" to
            specify mean or summed intensities, respectively.
        min_rel_intensity (float): minimum relative intensity for peaks to
            be considered. NOTE: must be decimal fraction in range 0 -> 1.
        min_abs_intensity (float): minimum absolute intensity for peaks to be
            considered.
        min_peak_identity (float): minimum % of spectra for peaks to be present
            in to be included in final consensus spectrum. NOTE: this must
            be expressed as a decimal fraction in range 0 -> 1.

    Returns:
        float: MA value for consensus spectrum (number of unique peaks after
            binning and intensity thresholding).
    """
    #  generate consensus spectrum
    consensus_spectrum = generate_consensus_spectrum(
        spectra=spectra,
        min_rel_intensity=min_rel_intensity,
        min_peak_identity=min_peak_identity,
        ppm_window=ppm_window,
        intensity_grouping=intensity_grouping,
        min_abs_intensity=min_abs_intensity
    )

    #  calculate and return MA of consensus spectrum
    return calculate_ma_spectrum_raw(spectrum=consensus_spectrum)
