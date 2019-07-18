"""Module for filtering spectra based on a selection of criteria
See spectra_filters.py for details

.. moduleauthor:: Graham Keenan 2019

"""

import json
from typing import Tuple, List, Union
import extractors.spectra_filter.spectra_filters as flt

# Constants
MS1, MS2 = "ms1", "ms2"
DEFAULT_ERROR = 0.01
PPM_ERROR = 1*(10**6)


class MissingKeysException(Exception):
    """Exception for missing MS1 and/or MS2 keys in ripper file
    """


def load_ripper_file(filename: str) -> Tuple[dict, dict]:
    """Loads a Ripper JSON file and obtians the MS1 and MS2 data

    Args:
        filename (str): Name of the file

    Raises:
        MissingKeysException: MS1 and/or MS2 key missing

    Returns:
        Tuple[dict, dict]: MS1 and MS2 data
    """

    with open(filename)as f_d:
        data = json.load(f_d)

    if MS1 in data.keys() and MS2 in data.keys():
        return data[MS1], data[MS2]

    raise MissingKeysException("Lacking MS1 and/or MS2 key!")


class SpectraFilterExtractor:
    """Class for filtering MS data based on a sleection of criteria

    Args:
        ripper_file (str): Name of the MS data file
        use_ppm_error (bool, optional): Use PPM error value. Defaults to False
        ppm_value (int, optional): PPM value to use for PPM error. Defaults to 10

    Raises:
        MissingKeysException: Missing MS1 and/or MS2 keys
    """

    def __init__(
            self,
            ripper_file: str,
            use_ppm_error: bool = False,
            ppm_value: int = 10,
        ):
        self.ppm_error = use_ppm_error
        self.ppm_value = ppm_value

        self.ms1, self.ms2 = load_ripper_file(ripper_file)


    def find_target(
            self,
            target: float,
            mass_list: List[float]
        ) -> Union[List[float], List]:
        """Finds a target match in a list of potential values within a given error
        within a given error
        Uses PPM error if use_ppm set, else uses default error

        Args:
            target (float): Target to search for
            mass_list (List[float]): List to search through

        Returns:
            Union[List[float], List]: Matches if any found, empty otherwise
        """

        # Set error based on self.ppm_error value
        error = self.ppm_value * PPM_ERROR if self.ppm_error else DEFAULT_ERROR

        # Get a list of masses that are within range of the target mass
        matches = filter(
            lambda x: x >= (target - error) and x <= (target + error),
            mass_list
        )

        # No matches, return empty list
        if not matches:
            {}

        # Return the matches
        return matches


    def filter_retention_time(
            self,
            ms_level: int,
            ret_time_range: List[float]
        ) -> dict:
        """Filters a spectra collection based on a retention time range

        Args:
            ms_level (int): MS level to filter
            ret_time_range (List[float]): Range of retention times to return

        Returns:
            dict: Filtered spectra or empty
        """

        if ms_level == 1:
            return flt.filter_retention_time(self.ms1, ret_time_range)
        if ms_level == 2:
            return flt.filter_retention_time(self.ms2, ret_time_range)

        print(f"Level: {ms_level} is not a supported MS level!")
        return {}


    def filter_parent_ion(self, parent_ions: List[float]) -> List[dict]:
        """Filters MS2 spectra based on a selection of parent ions

        Args:
            parent_ions (List[float]): Parent ions to search for

        Returns:
            dict: Spectra that contain the parent
        """

        return flt.filter_parent_ion(self.ms2, parent_ions, self.find_target)


    def filter_signature_ions(
            self,
            ms_level: int,
            sig_ions: List[float]
        ) -> dict:
        """Filters spectra based on a selection of signature ions

        Args:
            ms_level (int): MS level to filter
            sig_ions (List[float]): Signature ions to search for

        Returns:
            dict: Spectra that contain a signature or empty
        """

        if ms_level == 1:
            return flt.filter_signature_ions(
                self.ms1, sig_ions, self.find_target
            )

        if ms_level == 2:
            return flt.filter_signature_ions(
                self.ms2, sig_ions, self.find_target
            )

        print(f"Level: {ms_level} is not a supported MS level!")
        return {}


    def filter_total_intensity(
            self,
            ms_level: int,
            thresh: int
        ) -> dict:
        """Filters spectra based on their total intensity
        Spectra must exceed the threshold to be valid

        Args:
            ms_level (int): MS level to filter
            thresh (int): Total intensity threshold

        Returns:
            dict: Spectra that exceed threshold, or empty
        """

        if ms_level == 1:
            return flt.filter_total_intensity(self.ms1, thresh)

        if ms_level == 2:
            return flt.filter_total_intensity(self.ms2, thresh)

        print(f"Level: {ms_level} is not a supported MS level!")
        return {}


    def filter_max_intensity(
            self,
            ms_level: int,
            thresh: int
        ) -> dict:
        """Filters spectra absed on their maximum intensity
        Spectra must exceed the threshold to be valid

        Args:
            ms_level (int): MS level to filter
            thresh (int): Maximum intensity threshold

        Returns:
            dict: Spectra that exceed the threshold, or empty
        """

        if ms_level == 1:
            return flt.filter_max_intensity(self.ms1, thresh)

        if ms_level == 2:
            return flt.filter_max_intensity(self.ms2, thresh)

        print(f"Level: {ms_level} is not a supported MS level!")
        return {}


    def filter_mass_difference(
            self,
            ms_level: int,
            mass_diffs: dict,
            comparisons: int
        ) -> dict:
        """Filters spectra based on the mass difference between peaks.
        Scans through a list of monomor mass differences and checks for peaks
        that match the difference

        Args:
            ms_level (int): MS level to filter
            mass_diffs (dict): Monomers and the associated mass differences
            comparisons (int): Total number of comparisons to perform

        Returns:
            dict: Spectra IDs with found monomers
        """

        if ms_level == 1:
            return flt.filter_mass_difference(
                self.ms1, mass_diffs, comparisons, self.find_target
            )

        if ms_level == 2:
            return flt.filter_mass_difference(
                self.ms2, mass_diffs, comparisons, self.find_target
            )

        print(f"Level: {ms_level} is mot a supported MS level!")
        return {}
