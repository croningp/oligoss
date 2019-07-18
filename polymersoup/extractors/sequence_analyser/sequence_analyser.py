"""Analyses Mass Spec data to search for sequencing information

Processes MS1 and MS2 data, looking for fragment data within
KEEP IT CLEAN, BARRY

.. moduleauthor:: Graham Keenan 2019
.. moduleauthor:: David Doran 2019

"""

import os
import sequence_helpers as hlp

# Constants
MAPI_THRESHOLD = 90

class SequenceAnalyser(object):
    """Class for parsing MS data by looking through sequencing data

    Arguments:
        filename {str} -- Name of the root file
        msdata {dict} -- Mass Spec data
        seq_data {dict} -- Sequencing data

    """

    def __init__(
            self,
            filename: str,
            msdata: dict,
            seq_data: dict,
            output_dir: str,
            bruker=True
        ):
        self.filename = filename
        self.ms1 = msdata["ms1"]
        self.ms2 = msdata["ms2"]
        self.sequences = seq_data
        self.output_dir = os.path.abspath(output_dir)
        self.bruker = bruker


    def process_ms1(self, sequence_info: dict) -> (float, list, list):
        """Processes MS1 data

        Searches for sequence masses in the MS1 data
        Creates a combined EIC of all matches

        Arguments:
            sequence_info {dict} -- [description]

        Returns:
            float -- Maximum intensity of the most intense EIC
            list -- List of combined EICs
            list -- List of conformed masses in the MS1
        """

        sequence_masses = sequence_info["MS1"]
        combined_EICs, confirmed_masses = [], []

        # Go through each MS1 spectrum
        for spectrum in self.ms1.values():
            EICs = []
            mass_list = spectrum["mass_list"]

            # Search for sequence masses in the MS1 & create EIC
            for target in sequence_masses:
                matches = hlp.find_target_matches(target, mass_list, bruker=self.bruker)
                eic = hlp.make_EIC(matches, spectrum)

                # Only keep those that have an intensity associated
                if eic[1] > 0:
                    EICs.append(eic)
                confirmed_masses.append(target)

            # Remove duplicate retention times and construct combined EIC
            EICs = hlp.combine_duplicate_retention_times(EICs)
            if EICs:
                combined_EICs.append(EICs)

        # Get maximum intensity and convert confirmed masses to Set
        if combined_EICs:
            maximum_intensity = hlp.maximum_intensity_in_eics(combined_EICs)
        else:
            maximum_intensity = 0

        confirmed_masses = sorted(list(set(confirmed_masses)))

        return maximum_intensity, combined_EICs, confirmed_masses


    def process_ms2(
            self,
            ms1_confirmed_masses: list,
            sequence_info: dict
        ) -> (list, list):
        """Processes all MS2 data
        Goes through all fragment and unique fragment data, searching for masses

        Arguments:
            ms1_confirmed_masses {list} -- Masses confirmed in MS1
            sequence_info {dict} -- Sequencing Information

        Returns:
            list -- MS2 confirmed fragments
            list -- MS2 combined EIC
        """

        # Don't proces sif we have no MS1 confirmed massses
        if not ms1_confirmed_masses:
            return

        confirmed_fragments = []
        fragments = sequence_info["fragments"]
        sequence_peaks = sequence_info["peak_list"]
        unique_fragments = sequence_info["unique_fragments"]

        # Go through each MS2 spectrum
        for spectrum in self.ms2.values():
            # Search for fragments in MS2
            found_fragments = self.process_fragment_data(
                ms1_confirmed_masses,
                fragments,
                sequence_peaks,
                spectrum
            )

            # Add to list if fragments found
            if found_fragments:
                print(f"Found fragments: {found_fragments}")
                confirmed_fragments.extend(found_fragments)

        # Convert to Set
        confirmed_fragments = list(set(confirmed_fragments))

        if confirmed_fragments:
            print(f"Confirmed fragments: {confirmed_fragments}")

        # Create combined EIC for unqiue fragments
        confirmed_fragments, ms2_combined_eic = self.process_unique_fragments(
            unique_fragments,
            confirmed_fragments,
            fragments
        )

        return confirmed_fragments, ms2_combined_eic


    def process_fragment_data(
            self,
            ms1_masses: list,
            fragments: dict,
            peaks: list,
            spectrum: dict
        ) -> (list, str):
        """Looks for fragment information in the MS2 data

        If data is found, those fragments are confirmed and added to a list

        Arguments:
            ms1_masses {list} -- Confirmed MS1 masses
            fragments {dict} -- Fragment Information
            peaks {list} -- Fragment peaks
            spectrum {dict} -- MS2 spectra information

        Returns:
            list, str -- COnfirmed fragments if any, empty string otherwise
        """

        parent = spectrum["parent"]
        mass_list = spectrum["mass_list"]

        # Don't process if parent isn't present in MS1 confirmed masses
        if not hlp.check_ms2_parent_in_ms1(
                parent, ms1_masses, bruker=self.bruker
            ):
            return ""

        # Don't process if MAPI thresholf not met
        if not hlp.check_ms2_mapi(
                peaks, spectrum,
                bruker=self.bruker, mapi_threshold=MAPI_THRESHOLD
        ):
            return ""

        # Check fragments are in MS2
        confirmed_fragments = []
        for fragment, fragment_masses in fragments.items():
            fragment_present = hlp.check_fragment_in_ms2(
                fragment,
                fragment_masses,
                mass_list,
                bruker=self.bruker
            )

            if fragment_present:
                confirmed_fragments.append(fragment_present)

        if confirmed_fragments:
            return confirmed_fragments

        return ""


    def process_unique_fragments(
            self,
            unique_fragments: list,
            confirmed_fragments: list,
            fragment_info: dict
        ) -> (list, list):
        """Looks for unique fragment information in the MS2, if any

        If unique fragments, creates a combined EIC of those

        Arguments:
            unique_fragments {list} -- Unique Fragment Labels
            confirmed_fragments {list} -- List of confirmed fragments
            fragment_info {dict} -- Fragment information
        Returns:
            list, list -- Confirmed fragments and MS2 combined EIC
        """

        # Don't process if no unique fragments
        if not unique_fragments:
            return confirmed_fragments, []

        # Get unique fragments that are confirmed
        unique_confirmed_fragments = [
            frag for frag in confirmed_fragments if frag in unique_fragments
        ]

        # Don't process if no confirmed unique fragments
        if not unique_confirmed_fragments:
            return confirmed_fragments, []

        # Construct combined EIC
        combined_EIC = hlp.get_combined_unique_EIC(
            unique_confirmed_fragments,
            fragment_info,
            self.ms2.values(),
            bruker=self.bruker
        )

        return confirmed_fragments, combined_EIC


    def write_info(self, ms1_threshold: int):
        """Writes out common information used during the current run
        Written out in JSON format

        Args:
            ms1_threshold (int): Threshold for MS1 to meet before processing MS2
        """

        info = {
            "mapi_threshold": MAPI_THRESHOLD,
            "use_bruker_error": self.bruker,
            "peak_match_error": hlp.BRUKER_ERROR if self.bruker else hlp.DEFAULT_ERROR,
            "ms1_threshold": ms1_threshold
        }

        hlp.write_json(
            info,
            os.path.join(self.output_dir, "info.json")
        )


    def parse_data(self, ms1_threshold=1000):
        """Processes each sequence, looking for relevant information
        in the MS1 and MS2 data.

        If MS1 total intensity is less than the MS1 threshold, don't process MS2

        Keyword Arguments:
            ms1_threshold {int} -- Threshold for MS1 intensities (default: {1000})
        """

        ms1_output, ms2_unique_output, ms2_output = {}, {}, {}
        self.write_info(ms1_threshold)

        # Go through each sequence from the sequencing data
        for sequence, seq_info in self.sequences.items():
            print(f"Processing Sequence: {sequence}")
            seq_output = {}
            uniques = seq_info["unique_fragments"]

            # Process MS1
            ms1_max, ms1_combined_eic, ms1_confirmed_targets = self.process_ms1(
                seq_info
            )

            # Don't process MS2 if most intense MS1 below threshold
            if ms1_max < ms1_threshold:
                continue

            # Add MS1 data to output
            ms1_output[sequence] = ms1_combined_eic

            # Process MS2
            confirmed_fragments, ms2_combined_eic = self.process_ms2(
                ms1_confirmed_targets,
                seq_info
            )

            # Check if we have confirmed fragments
            if confirmed_fragments:
                print(f"MS2 -- Confirmed Fragments: {confirmed_fragments}")

                # Add confirmed to output
                seq_output["confirmed_fragments"] = confirmed_fragments

                # Use MS2 EIC if unqiue fragments, else use MS1 EIC
                seq_output["EIC"] = ms2_combined_eic if uniques else ms1_combined_eic

                # Get confirmed unique fragments
                unique_in_confirmed = [x for x in uniques if x in confirmed_fragments]

                if uniques and unique_in_confirmed:
                    # Add to unique output if we have uniques and are confirmed
                    ms2_unique_output[sequence] = seq_output
                else:
                    # Usual MS2 output
                    ms2_output[sequence] = seq_output

        # Write out to file
        hlp.write_json(
            ms1_output,
            os.path.join(self.output_dir, f"{self.filename}_ms1_eics.json")
        )

        hlp.write_json(
            ms2_output,
            os.path.join(self.output_dir, f"{self.filename}_ms2_info.json")
        )

        hlp.write_json(
            ms2_unique_output,
            os.path.join(self.output_dir, f"{self.filename}_ms2_unique_info.json")
        )
