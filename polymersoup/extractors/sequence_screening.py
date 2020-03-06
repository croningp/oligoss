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
    (ripper_dict).

    Args:
        ions (list of floats): list of target ion m/z values.
        ms_level (int): MS level at which to screen for ions.
        ripper_dict (dict): mzml ripper mass spectrometry data dictionary.
        err (float): error tolerance for matching target ions to masses found
                in spectra; units either in absolute mass units or ppm.
        err_abs (bool, optional): specifies whether err is in absolute mass
                units or ppm; if True, err is in absolute mass units.
                Defaults to True.
        min_max_intensity (float, optional): minimum maximum intensity of EIC
                for it to be returned. Defaults to None.
        min_total_intensity (float, optional): maximum total intensity of EIC
                for it to be returned. Defaults to None.
        precursors (list of floats, optional): list of precursor ions for
                matching to parents/precursors of MSn spectra.
        mapi (optional, float): minimum annotated peak intensity (in % of most
                intense peak), specifies minimum relative intensity of most
                intense peak associated with a sequence for spectra to be used
                in screening.
        peak_list (optional, list of floats): list of ALL ions associated with
                targets, not just ions used in generating EIC.
    Returns:
        EIC (list of lists): extracted ion chromatogram in format: [[Rt, I]..]
                where Rt = retention time (float) and I = intensity of all
                species that match target ions at that retention time.
    """
    ripper_dict_at_ms_level = ripper_dict[f'ms{ms_level}']
    eic = []

    # Get error thresholds
    if err_abs:
        errs = [float(err) for ion in ions]
    else:
        errs = [((ion * err) * 10**-6) for ion in ions]

    # Establish target thresholds before entering loops
    target_ion_limits = [
        (ion, ion - errs[i], ion + errs[i])
        for i, ion in enumerate(ions)
    ]

    # Sum intensities
    for spectrum_info in ripper_dict_at_ms_level.values():
        intensity = 0
        for mass in spectrum_info['mass_list']:
            for ion, min_target, max_target in target_ion_limits:
                if min_target <= mass <= max_target:
                    mass_key = f'{mass:.4f}'
                    intensity += spectrum_info[mass_key]
        if intensity > 0:
            eic.append(
                [float(spectrum_info['retention_time']), float(intensity)]
            )

    eic = sorted(eic, key=lambda item: item[1])

    # Return empty list if max intensity is below min_max_intensity
    if min_max_intensity and eic:
        # EIC already sorted so just check last element
        if eic[-1][1] < min_max_intensity:
            return []

    # Return empty list if sum of intensities is less than min_total_intensity
    if min_total_intensity and eic:
        if sum([item[1] for item in eic]) < min_total_intensity:
            return []

    return eic


def generate_EIC_deprecated(
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
    (ripper_dict).

    Args:
        ions (list of floats): list of target ion m/z values.
        ms_level (int): MS level at which to screen for ions.
        ripper_dict (dict): mzml ripper mass spectrometry data dictionary.
        err (float): error tolerance for matching target ions to masses found
                in spectra; units either in absolute mass units or ppm.
        err_abs (bool, optional): specifies whether err is in absolute mass
                units or ppm; if True, err is in absolute mass units.
                Defaults to True.
        min_max_intensity (float, optional): minimum maximum intensity of EIC
                for it to be returned. Defaults to None.
        min_total_intensity (float, optional): maximum total intensity of EIC
                for it to be returned. Defaults to None.
        precursors (list of floats, optional): list of precursor ions for
                matching to parents/precursors of MSn spectra.
        mapi (optional, float): minimum annotated peak intensity (in % of most
                intense peak), specifies minimum relative intensity of most
                intense peak associated with a sequence for spectra to be used
                in screening.
        peak_list (optional, list of floats): list of ALL ions associated with
                targets, not just ions used in generating EIC.
    Returns:
        EIC (list of lists): extracted ion chromatogram in format: [[Rt, I]..]
                where Rt = retention time (float) and I = intensity of all
                species that match target ions at that retention time.
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
    MS1 EICs for sequences.

    Args:
        silico_dict (dict): full MSMS in silico sequence dict.
        ripper_dict (dict): mzml ripper mass spectrometry data dictionary.
        err (float): error tolerance for matching target ions to masses found
                in spectra; units either in absolute mass units or ppm.
        err_abs (bool, optional): specifies whether err is in absolute mass
                units or ppm; if True, err is in absolute mass units.
                Defaults to True.
        min_max_intensity (float, optional): minimum maximum intensity of EIC
                for it to be returned. Defaults to None.
        min_total_intensity (float, optional): maximum total intensity of EIC
                for it to be returned. Defaults to None.
    """
    for sequence in silico_dict:
        if type(silico_dict[sequence]) == dict:
            masses = silico_dict[sequence]["MS1"]
        else:
            masses = silico_dict[sequence]

        silico_dict[sequence] = [float(mass) for mass in masses]

        sequence_MS1_EIC = generate_EIC(
            ions=masses,
            ms_level=1,
            ripper_dict=ripper_dict,
            err=err,
            err_abs=err_abs,
            min_max_intensity=min_max_intensity,
            min_total_intensity=min_total_intensity
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
    returns list of confirmed fragments associated with precursors.

    Args:
        precursors (list of floats): list of m/z values for sequence precursors.
        fragment_dict (dict): dictionary of fragments and associated masses
            in insilico fragment dict format.
        peak_list (list of floats): list of ALL MS1 and MS2 / MSn masses
            associated with sequence.
        ripper_dict (dict): dict of mass spec data in mzml ripper format.
        err (float): error tolerance used when matching masses; units can
            either be in absolute mass units or ppm.
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
            mass spectrum data from ripper_dict.
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
    m/z value of most intense peak and I = intensity of most intense peak.

    Args:
        spectrum (dict): single spectrum dictionary in mzml ripper format.

    Returns:
        most_intense (tuple): tuple of peak m/z and intensity in format
                    (m/z, I).
    """

    masses = []
    
    # if there are no masses in spectrum, return 0 
    if not spectrum['mass_list']:
        return (0,0)

    # get list of mass (m/z), intensity tuples for all peaks in spectrum
    for mass in spectrum['mass_list']:
        mass_key = str(mass) + "0000"
        mass_key = mass_key[0:mass_key.find('.')+5]
        masses.append(
            (float(mass), float(spectrum[mass_key]))
        )

    # sort peaks by most intense
    masses = sorted(masses, key = lambda x: x[1])

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
    intensity of the most intense matching peak.

    Args:
        spectrum (dict): spectrum dict in mzml ripper dict format.
        peak_list (list): list of theoretical m/z values associated with a
                        sequence.
        err (float): error tolerance used when matching masses; units can
            either be in absolute mass units or ppm.
        err_abs (bool, optional): specifies whether err units are in absolute
            mass units or ppm; if True, absolute mass unit is used.
            Defaults to True.

    Returns:
        top_intensity (float): intensity of most intense matching peak.
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
    to most intense peak in the spectrum.

    Args:
        spectrum (dict): spectrum dict in mzml ripper format
        peak_list (list of floats): list of peaks (m/z values), typically peaks
                    associated with a target sequence.
        err (float): error tolerance in matching peak masses - can be in
                    absolute mass units or ppm.
        err_abs (bool): specifies whether err units are in absolute mass units
                    or ppm; if True, units = absolute mass units.

    Returns:
        float: maximum annotated peak intensity in % of most intense peak.
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
    fragments that have matching ions in the spectrum.

    Args:
        spectrum (dict): spectrum dict in mzml ripper format.
        fragment_dict (dict): fragment dictionary in standard insilico format.
        err (float): error tolerance for matching fragment masses; units can
                    either be absolute mass units or ppm.
        err_abs (bool): specifies whether err units are in absolute mass
                    units or ppm; if True, err is in absolute mass units.
        essential_signtures (list): list of signature fragment ids that must
                    be found in spectrum if any fragments are to be confirmed
                    from said spectrum (defaults to []).

    Returns:
        confirmed_fragments (list of strings): list of fragment ids from
                    input fragment_dict that have one or more matching masses
                    in the spectrum.
    """
    # initiate list to store confirmed fragments
    confirmed_fragments = []

    # retrieve masses from spectrum
    spectrum_masses = [float(mass) for mass in spectrum['mass_list']]

    # iterate through SIGNATURE FRAGMENTS
    for signature, signature_masses in fragment_dict['signatures'].items():
        
        if type(signature_masses) == list:
           
            # retrieve masses in spectrum that match fragment target masses
            matches = find_multiple_targets(
                targets=signature_masses,
                candidates=spectrum_masses,
                err=err,
                err_abs=err_abs
            )

            # if any spectrum masses are a match for fragment masses, add
            # fragment id to list of confirmed fragments
            if matches:
                confirmed_fragments.append(signature)

    # check for terminal modification signatures
    if fragment_dict['signatures']['terminal_modifications']:
        
        # iterate through terminal modifications to find free modification 
        # fragments
        for mod, mod_ions in fragment_dict['signatures']['terminal_modifications'].items():

            matches = find_multiple_targets(
                targets=mod_ions,
                candidates=spectrum_masses,
                err=err,
                err_abs=err_abs
            )
            if matches:
                confirmed_fragments.append(mod)

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
                targets=masses,
                candidates=spectrum_masses,
                err=err,
                err_abs=err_abs
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
    can be confirmed as associated with those sequences from the ripper data.

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
        associated with sequence, including signature ions.

        ripper_dict (dict): mass spec data dict in mzml ripper format.
        err (float): error tolerance for screening masses.
        err_abs (bool): specifies whether units of err are in absolute mass
        units or ppm; if True, units = absolute mass units.
        min_ms2_peak_abundance (float): minimum relative abundance (in % of
        most abundant / intense MS2 spectrum peak) of most intense peak in
        peak_list.
        essential_signatures (list of strings): list of signatures that
        MUST be found in a spectrum for fragments to be confirmed from that
        spectrum.
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

def generate_bpc(
    msdata,
    N_peaks_per_spectrum,
    ms_level=1,
    min_intensity=0,
    dymanic_exclusion_window=None,
    abs_precursor_exclusion=False
):
    """
    Screens mass spectra dict in mzml_ripper format and returns base peak 
    chromatogram (bpc) at specified ms level, along with info on most intense
    peaks in data
    
    Args:
        msdata (dict): dict of mass spectra in mzml ripper format. NOTE: this
            can either be dict of ALL ms levels or specified target ms level
        N_peaks_per_spectrum (int): number of most intense peaks to keep 
            track of per spectrum for  'intense_peaks' list in bpc_dict. 
        ms_level (int, optional): specifies which ms level is to be used to
            generate bpc. Defaults to 1.
        min_intensity (int, optional): minimum absolute intensity for peaks to 
            be added to 'intense_peaks' list in bpc dict. Defaults to 0.
    
    Returns:
        dict: bpc dict in format: 
            {   
                'intense_peaks': [[mz, I]...]
                'spectrum_id': {
                    'retention_time': retention_time,
                    'peaks' [[mz, I]...]
                },
                ...
            }
            where mz, I = m/z, intensity of peaks. 
            N most intense peaks of each spectrum (N is specified by 
            N_peaks_per_spectrum arg) are stored for each spectrum and added
            to both spectrum subdict and 'intense_peaks' list 
    """

    # remove all data in ripper dict except relevant MS level
    if f'ms{ms_level}' in msdata:
        msdata = msdata[f'ms{ms_level}']

    # init list to store most intense peaks from each spectrum
    intense_peaks = []

    # init dict to store bpc
    bpc_dict = {}

    # iterate through spectra, retrieving most intense peaks
    for spectrum_id, spectrum_info in msdata.items():

        # check whether current dict item is a spectrum with relevant info
        if not spectrum_id.startswith('spectrum_'):
            pass
        else: 
            
            # find N most intense peaks in current spectrum and add to list
            # of most intense peaks in format [[mz, I]...] where mz, I = 
            # m/z and intensity of peaks, respectively 
            spectrum_peaks = return_most_intense_peaks_spectrum(
                spectrum=spectrum_info,
                N_peaks=N_peaks_per_spectrum)

            # add spectrum peaks to most intense peaks 
            intense_peaks.extend(spectrum_peaks)

            # remove duplicate peaks and ensure all peaks meet minimum 
            # intensity threshold
            intense_peaks = list(set(
                [peak for peak in intense_peaks
                if float(peak[1]) >= min_intensity]))
            
            # add spectrum retention time and peaks to bpc_dict
            bpc_dict[spectrum_id] = {
                'retention_time': spectrum_info['retention_time'],
                'peaks': spectrum_peaks
            }
    
    # add most intense peaks across all data to bpc
    bpc_dict['intense_peaks'] = intense_peaks
    
    return bpc_dict

def fingerprint_screen_MSn_from_bpc_precursors(
    bpc_dict,
    msdata,
    silico_info,
    err,
    err_abs,
    ms_level=2,
    comparisons_per_spectrum=None
):
    """
    This function takes info from a base peak chromatogram (bpc) and screens 
    MSn spectra for fragments of precursors in bpc that match fingerprints for
    specific monomers. Fingerprints can be for mass differences between MSn 
    product ion peaks that correspond to potential monomer mass ladders or
    single MSn product ion peaks that correspond to monomer-specific signature 
    ions. 

    Args:
        bpc_dict (dict): base peak chromatogram info. For format see function: 
            'generate_bpc' (in sequence_screening.py)
        msdata (dict): mass spectra in mzml ripper format 
        silico_info (dict): dictionary of monomers, associateed massdiffs and 
            signature ions 
        err (float): error threshold for matching signature ions and massdiffs 
            in spectra. Units can either be in ppm or absolute mass units 
        err_abs (bool): specifies whether err is in units of ppm or absolute
            mass units; if True, units = absolute mass units; if False, units =
            ppm
        ms_level (int, optional): specifies ms level to screen for massdiffs
            and signatures in product ions. NOTE: this MUST be greater than 1.
            Defaults to 2. 
        comparisons_per_spectrum (int, optional): specifies N most intense peaks
            to use for massdiff comparison in target spectra. Defaults to 20.
    
    Returns:
        dict: dict of spectral ids and associated hits for monomer signature 
            ions and / or massdiffs. Format: 
            {
                'spectrum_id': {
                    'X': {
                        'mass_diffs': [m/z,...],
                        'signatures': [m/z,...]
                    }
                },
                ...
            }
        where 'X' = monomer one-letter code, m/z = m/z value of peaks that 
        match signatures for 'X' or are separated from other peaks by massdiff
        of 'X'
    """

    # check whether bpc info has been input; if so, use bpc to select precursors
    # for screening. Otherwise, screen every precursor in MSn spectra 
    if bpc_dict:

        # get list of precrusor ions of interest for filtering MSn spectra further 
        precursors = [
            float(peak[0]) for peak in bpc_dict['intense_peaks']]
    
        # remove MSn spectra that do not have precursors of interest 
        filtered_msdata = filter_parent_ion(
            msdata=msdata,
            parent_ions=precursors,
            err=err,
            err_abs=err_abs,
            ms_level=ms_level)
    
    else:
        filtered_msdata = {
            ms: msinfo for ms, msinfo in msdata.items()
            if ms == f'ms{ms_level}'}
    
    # init dict to store spectral hits for monomer fingerprints, and add silico
    # info containing theoretical monomer signatures and massdiffs for reference]
    monomer_hits = {'silico_info': silico_info}

    # screen MSn spectra for monomer fingerprints and / or massdiffs 
    monomer_hits.update(filter_monomer_fingerprints(
        msdata=filtered_msdata[f'ms{ms_level}'],
        monomer_massdiffs=silico_info,
        total_comparisons=comparisons_per_spectrum,
        err=err,
        err_abs=err_abs,
        ms_level=ms_level))
    
    # return monomer massdiff / signature ion hits 
    return monomer_hits

def screen_MS1_precursors_masslib(
    MS1_masslib,
    extractor_parameters,
    spectral_assignments
):
    print(f'extractor keys = {extractor_parameters.keys()}')
    print(f'extractor_parameters error = {extractor_parameters["error"]}')
    err = float(extractor_parameters["error"])
    err_abs = extractor_parameters["err_abs"]

    spectral_assignments = {
        spectrum: spectrum_info
        for spectrum, spectrum_info in spectral_assignments.items()
        if spectrum.find('silico') == -1
    }

    sequence_precursor_assignments, assigned_precursors = {}, []

    for sequence, masses in MS1_masslib.items():

        for spectrum, spectral_assignment in spectral_assignments.items():

            precursor = spectral_assignment["precursor"]
            
            seq_match = find_target(
                target=precursor,
                candidates=masses,
                err=err,
                err_abs=err_abs)
            
            if seq_match: 
                
                assigned_precursors.append([spectrum, precursor])
                if sequence not in sequence_precursor_assignments:

                    sequence_precursor_assignments[sequence] = [
                        [spectrum, precursor]]
                
                elif sequence in sequence_precursor_assignments: 

                    sequence_precursor_assignments[sequence].append([spectrum, precursor])
    if sequence_precursor_assignments:
        sequence_precursor_assignments.update({
            'assigned_precursors': assigned_precursors,
            'MS1_silico_lib': {
                seq: MS1_masslib[seq]
                for seq in MS1_masslib
                if seq in sequence_precursor_assignments
            }
        })
    
    
    return sequence_precursor_assignments 

def extract_MS1(
    silico_dict,
    extractor_parameters,
    ripper_dict,
    return_Rt_list=True
):
    """
    Takes an in silico sequence dict, dictionary of extractor parameters
    (passed on from input parameters .json file), mass spec data in mzml ripper
    format, and outputs a dictionary of sequences and their corresponding
    MS1 EICs (if there are MS1 hits).

    Args:
        silico_dict (dict): dictionary of sequences and corresponding masses.
        extractor_parameters (dict): dictionary of extractor parameters in
            standard input parameters format.
        ripper_dict (dict): dictionary of mass spec data in mzml ripper format.
        return_Rt_list (bool, optional): specify whether to return a list of 
            ALL retention times in the ripper dict. Defaults to True. 

    Returns:
        EIC_dict (dict): dictionary of sequences and corresponding MS1 EICs.
    """
    # check whether file path to input parameters has accidentally been passed
    # in instead of extractor_parameters dict
    if type(extractor_parameters) == str:
        extractor_parameters = read_parameters(
            extractor_parameters)['extractor_parameters']

    # check whether full input parameters dict has accidentally been passed 
    # in instead of extractor parameters dict 
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

    # write list of all retention times in ripper data to output dict 
    if return_Rt_list: 
        
        retention_times = []

        for spectrum_info in ripper_dict['ms1'].values():
            retention_times.append(float(spectrum_info['retention_time']))

        EIC_dict['retention_time'] = sorted(retention_times)

    # return dict of found sequences and their corresponding MS1 EICs
    return EIC_dict

def extract_MS2(
    silico_dict,
    extractor_parameters,
    ripper_dict,
    filter_most_intense=True,
    extract_uniques=True,
    return_Rt_list=True
):
    """
    This function takes a silico MSMS sequence dict and generates MS2 EICs
    for each sequence in the dictionary, returning a dictionary of sequences
    and corresponidng MS2 EICs.

    Args:
        silico_dict (dict): full MSMS silico dict in standard format.
        extractor_parameters (dict): extractor parameters dict passed on
            from input parameters .json file.
        ripper_dict (dict): mass spec data dict in mzml ripper format.
        filter_most_intense (bool, optional): specifies whether to only
            return the most intense MS2 EIC for target sequences.
            Defaults to True.

    Returns:
        EIC_dict: dictionary of sequences and their corresponding MS2 EICs.
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

def find_monomer_fingerprints_MSn_spectrum(
    spectrum,
    monomer_fingerprints, 
    min_ms2_intensity,
    err,
    err_abs,
    N_peak_comparisons=None
):

    # get spectrum peaks, sorted in descending order by peak intensity (i.e.
    # most intense first)
    peaks = return_most_intense_peaks_spectrum(
        spectrum=spectrum,
        N_peaks=N_peak_comparisons,
        min_relative_intensity=min_ms2_intensity)
    
    # get base peak / dominant ion (i.e. most intense peak in spectrum)
    most_intense_peak = peaks[0]
    
    # init dict to store fingerprint massdiff and / or signature for target
    # monomers 
    monomer_assignments = {}
    
    # set basepeak assignment to False. This will change to True if basepeak/
    # dominant ion is found in one or more monomer signature and / or mass 
    # ladder 
    basepeak_assignment = False

    # iterate through monomers and their fingerprints, and identify any 
    # found in spectrum 
    for monomer, fingerprint_info in monomer_fingerprints.items():
        
        # init list to store signature matches, and populate if any are found
        signature_matches = []
        
        if fingerprint_info['signatures']:
            signature_matches = find_multiple_targets(
                    targets=fingerprint_info['signatures'],
                    candidates=[peak[0] for peak in peaks],
                    err=err,
                    err_abs=err_abs)
        
        # get massdiffs for monomer and check for these in spectrum 
        massdiffs_monomer = fingerprint_info['massdiffs']
        massdiff_hits = find_massdiffs_spectrum_peaks(
            spectrum_peaks=peaks,
            massdiffs=massdiffs_monomer,
            err=err,
            err_abs=err_abs)
        
        # if matches for signatures and / or massdiffs have been found in 
        # spectrum, add spectrum peak assignments for monomer to 
        # monomer_assignments dict 
        if signature_matches or massdiff_hits: 
            
            massdiff_hits = [round(hit, 4) for hit in massdiff_hits]
            signature_matches = [round(match, 4) for match in signature_matches]

            # check whether most intense ion / basepeak has been assigned to
            # monomer mass ladder(s) and / or signature ions
            if round(float(most_intense_peak[0]), 4) in massdiff_hits or signature_matches:
                basepeak_assignment=True
           
            # get full assignment info for monomer fingerprints in spectrum, 
            # including basepeak assignments 
            monomer_assignments[monomer] = {
                'signatures': signature_matches,
                'mass_diffs': massdiff_hits
            }

    if monomer_assignments:
        monomer_assignments.update(
            {
                'base_peak': most_intense_peak,
                'base_peak_assignment': basepeak_assignment
            }
        )
    # return monomer assignments 
    return monomer_assignments 

def filter_monomer_fingerprints(
        msdata: dict,
        monomer_massdiffs: dict,
        total_comparisons: int,
        err: float,
        err_abs=True,
        ms_level=2,
        single_spectrum=False,
        single_spectra_id=None,
        find_common_peaks=True,
        min_ms2_intensity=0.1
    ) -> dict:
    """ Filters the spectra based on the mass difference between peaks.
    Scans through a list of monomer mass differences and check for peaks that
    match the difference.

    Args:
        msdata (dict): full mass spec data JSON (in mzml ripper format) OR
            subset of mzml ripper dict for one ms level.
        monomer_massdiffs (dict): monomers and their associated mass differences.
        total_comparisons (int): total number of comparisons to perform.
        target_func (Callable): find target function.
        err (float): error threshold for screening - either in absolute mass
            units or ppm.
        err_abs (bool): specifies whether err units are absolute mass units
            or ppm.
        ms_level(int):specifies ms_level for spectra being screened (default: 2).
        single_spectrum (bool): specifies if data is single spectrum or not
            (default: False).
        spectra_id (str): spectra id string of single spectrum data (default: None).
        find_common_peaks (bool, optional): specifies whether to identify peaks
            that are common to two or more mass ladders for monomers (default: 
            True).
        min_ms2_intensity (float, optional): minimum MS2 peak intensity as a 
            decimal fraction of dominant ion. Peaks with an intensity below
            this value will be discarded from consideration. (default: 0.01).
            NOTE: must be in range 0 -> 1. 

    Returns:
        dict: Spectra IDs with found monomers.
    """

    # checks whether input is full ripper dict or just one ms level
    if f'ms{ms_level}' in msdata:
        msdata = msdata[f'ms{ms_level}']

    # init dict to store spectral assignments for monomer fingerprints 
    assigned_spectra = {}

    # get info on monomer fingerprints, excluding list of expected dominant
    # signatures 
    monomer_fingerprint_info = {
            k: v 
            for k,v in monomer_massdiffs.items()
            if k.find('dominant') == -1}

    # Iterate through each spectrum
    for spec_id, spectrum in msdata.items():
        
        monomer_fingerprint_hits = find_monomer_fingerprints_MSn_spectrum(
            spectrum=spectrum,
            monomer_fingerprints=monomer_fingerprint_info,
            err=err,
            err_abs=err_abs,
            min_ms2_intensity=min_ms2_intensity,
            N_peak_comparisons=total_comparisons)
        
        
        if monomer_fingerprint_hits:
            monomer_fingerprint_hits.update({'precursor': spectrum['parent']})
            assigned_spectra[spec_id] = monomer_fingerprint_hits
    # check whether common peak searching is specified; if so, update spectral
    # assignments with lists of peaks that are common to two or more mass 
    # ladders for target monomers
    if find_common_peaks: 
        return find_common_peaks_massdiff_spectra(
            spectral_assignments=assigned_spectra)

    # Return results
    return assigned_spectra