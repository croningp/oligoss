"""
This file contains classes for storing information on polymer-specific
features that are passed in from parameter handlers and passed on to
silico modules
"""
from ...utils.parameter_handlers.polymer_param_handlers import load_polymer_info
from ...utils.global_chemical_constants import FUNCTIONAL_GROUPS
from ...utils.errors import (
    InvalidMonomerId,
    InvalidMSFragmentation,
    InvalidModificationTarget
)
from ...utils.parameter_handlers.instrument_handlers import (
    retrieve_mass_spec_info_params_obj
)

class Polymer():
    def __init__(
        self,
        params_obj
    ):
        """
        Polymer object for use in silico functions.

        Args:
            params_obj (Parameters): instance of Parameters class

        Props:
            monomers (List[Monomer]): list of Monomer objects.
            dominant_signatures (List[str]): list of monomer id strings for
                monomers with dominant signature ions.
            mass_diff (float): mass lost upon addition of monomer to elongating
                backbone.
            elongation_unit (int): number of monomers added at each elongation.
            symmetry (bool): specifies whether Polymer termini are equivalent.
                True if equivalent, False if not.
            chain_terminators (): N/A.
            fragment_info (Dict[str, dict]): dict of fragment string ids and
                associated info.
            modifications (Dict[str, List[str]]): dict of sidechains and termini
                target for covalent modification(s).
        """

        #  set internal properties, params and polymer_config
        self.__params = params_obj
        self.__polymer_config = load_polymer_info(params_obj.polymer_class)

        #  set remaining properties
        self.monomers = self.instantiate_monomers()
        self.dominant_signatures = self.retrieve_dominant_signatures()
        self.mass_diff = self.retrieve_mass_diff()
        self.elongation_unit = self.__polymer_config["ELONGATION"]
        self.symmetry = self.__polymer_config["SYMMETRY"]
        self.chain_terminators = self.__polymer_config["CHAIN_TERMINATORS"]
        self.fragment_info = self.retrieve_fragment_info()

        self.modifications = self.build_modifications_dict(
            targets=self.__params.silico.modifications)

    def instantiate_monomers(self):
        """
        Returns monomers parameter for instance of Polymer class

        Returns:
            List[Monomer]: list of Monomer objects
        """
        monomer_ids = self.__params.monomers
        monomers = [
            Monomer(
                monomer_id=monomer_id,
                polymer_info=self.__polymer_config,
                mode=self.__params.mode,
                ms2_signature_types=self.__params.silico.ms2.signatures)
            for monomer_id in list(set(monomer_ids))]
        return monomers

    def retrieve_mass_diff(self):

        mass_diff = self.__polymer_config["MASS_DIFF"]
        try:
            mass_diff = str(mass_diff)
        except ValueError:
            if mass_diff[0] == "-":
                mass_diff = -FUNCTIONAL_GROUPS[mass_diff[1::]]
            else:
                mass_diff = FUNCTIONAL_GROUPS[mass_diff]

        return mass_diff

    def retrieve_fragment_info(self):
        """
        retrieves fragment_info dict from polymer config

        Returns:
            Dict[str, dict]: fragment_info dict
        """

        #  get list of fragments from input parameters
        fragment_series = self.__params.silico.ms2.fragment_series

        #  retrieve info for each fragment from polymer_config
        fragment_info = {
            fragment: self.__polymer_config["FRAG_SERIES"][fragment]
            for fragment in fragment_series
        }

        #  convert fragment mass_diffs to floats
        for series in fragment_info:
            massdiff_info = fragment_info[series]["mass_diff"]
            for mode, value in massdiff_info.items():
                if value and mode != "exceptions":
                    try:
                        fragment_info[series][mode] = float(value)
                    except ValueError:
                        if value[0] == "-":
                            fragment_info[series]["mass_diff"][mode] = (
                                -FUNCTIONAL_GROUPS[value[1::]]
                            )
                        else:
                            fragment_info[series]["mass_diff"][mode] = (
                                FUNCTIONAL_GROUPS[value])

        #  ensure all fragment series are compatible with instruments
        if not self.check_fragment_instrument_compatibility(fragment_info):
            raise InvalidMSFragmentation(
                frag_method=retrieve_mass_spec_info_params_obj(
                    params_obj=self.__params))

        return fragment_info

    def check_fragment_instrument_compatibility(self, fragment_info):
        """
        Checks whether proposed fragment series are compatible with the
        mass spec fragmentation methods available in instrument.

        Args:
            fragment_info (dict): dict of fragment series and associated
                info.

        Returns:
            bool: False if fragment series are not compatible with instrument,
                True if there is nothing to suggest fragment series and
                instrument are incompatible.
        """

        #  if no instrument has been specified, cannot check for incompatibility
        if not self.__params.instrument:
            return True

        #  retrieve MS2 fragmentation methods for instrument
        fragmentation_methods = retrieve_mass_spec_info_params_obj(
            params_obj=self.__params)["fragmentation"]["ms2"]

        #  init list to store permissible fragment series for instrument
        permissible_fragments = []

        #  iterate through fragmentation methods available in instrument,
        #  check if they are in the polymer_config
        for frag_method in fragmentation_methods:

            #  no default fragments are specified in polymer_config, so cannot
            #  check for incompatibility
            if "default_linear" not in self.__polymer_config["FRAG_SERIES"]:
                return True

            #  update permissible fragments if frag_method has compatible
            #  fragments for polymer class
            if frag_method in self.__polymer_config["FRAG_SERIES"][
                    "default_linear"]:
                permissible_fragments.extend(
                    self.__polymer_config["FRAG_SERIES"][
                        "default_linear"][frag_method])

        #  check that each fragment proposed is compatible with instrument
        #  fragmentation methods
        for frag in fragment_info:
            if frag not in permissible_fragments:
                raise Exception(frag, permissible_fragments)

        return True

    def retrieve_dominant_signatures(self):
        """
        Get list of monomer str codes for monomers expected to have dominant
        ms2 signatures

        Returns:
            List[str]: list of monomer str codes
        """
        dominant_sigs = []
        if not self.__params.silico.ms2.signatures:
            return None

        for sig in self.__params.silico.ms2.signatures:
            for monomer in self.monomers:
                if monomer.id in self.__polymer_config[
                        "MS_SIGNATURE_IONS"][sig]["dominant"]:
                    dominant_sigs.append(monomer.id)
        return dominant_sigs

    def build_modifications_dict(self, targets):
        """
        Takes targets of covalent modifications and builds dict of monomer ids
        and associated Modification objects

        Args:
            targets (Dict[str, List[str]]): dict of targets (either monomer_ids
                or termini ids) and lists of modification string ids that target
                monomer sidechains or termini.

        Returns:
            Dict[str, List[Modification]]: dict of monomer_ids or termini_ids
                and associated lists of Modification objects.
        """

        #  no targets so return empty dict
        if not targets:
            return {}

        #  build mod dict and return
        mod_dict = {
            str(target): [
                Modification(
                    polymer_info=self.__polymer_config,
                    mod_id=mod,
                    mode=self.__params.mode
                )
                for mod in targets[target]
            ]
            for target in targets
        }

        #  make sure all modifications have valid target
        for target, mods in mod_dict.items():
            if mods:
                for mod in mods:
                    valid_targets = []
                    if mod.side_chain_attachments:
                        valid_targets += mod.side_chain_attachments
                    if valid_targets and mod.termini:
                        valid_targets += [str(x) for x in mod.termini]
                    if valid_targets and str(target) not in valid_targets:
                        raise InvalidModificationTarget(
                            modification=mod,
                            target=target
                        )

        return mod_dict


class Monomer():
    def __init__(
        self,
        monomer_id,
        polymer_info,
        mode,
        ms2_signature_types
    ):
        """
        Defines Monomer objects to be used in generation of oligomer sequences.

        Args:
            monomer_id (str): id string of monomer (typically one letter code).
            polymer_info (dict): polymer-specific config dict.
            mode (str): denotes mass spec mode - "pos" or "neg" for positive and
                negative, respectively.
            ms2_signature_types (List[str]): list of signature types to
                determine which signature fragments are associated with Monomer.

        Properties:
            id (str): monomer id string (typically one letter code).
            neutral_mass (float): neutral monoisotopic mass of monomer.
            func_groups (List[list]): list containing monomer functional groups
                in format (using a standard amino acid as an example):
                    [["func_group_1", x], ["func_group_2", y]] where x, y =
                number of func_group_1, func_group_2 in monomer, respectively
                Amino acid example:
                    [["amine", 1], ["carboxyl", 1]].
            full_name (str): full name of monomer.
            loss_products (List[float or str]): list of neutral losses
                associated with monomer, denoted either by functional string
                or mass.
            ionizable (bool): True if monomer sidechain can ionize (i.e. take
                on an additional charge other than backbone charge), False if
                not.
            ion_info (List[str, float, float], optional): specifies exchangeable
                ion and charge range for ionizable sidechain. Format:
                    [{ion}, x, y] where {ion}, x, y = exchangeable ion, minimum
                        charge, maximum charge.
            signatures (List[float], optional): list of associated signature
                ion m/z values for monomer, if applicable.
        """

        #  set internal attr for polymer config
        self.__polymer_config = load_polymer_info(polymer_info)

        #  set props, inclusing "basic attributes" (neutral_mass, func_groups,
        #  full_name)
        self.__mode = mode
        self.id = monomer_id
        self.set_basic_attrs()

        # set remaining props
        self.loss_products = self.retrieve_loss_products()
        self.ionizable, self.ion_info = self.retrieve_sidechain_ions()
        self.signatures = self.retrieve_ms2_signature_ions(
            signature_types=ms2_signature_types)

    def set_basic_attrs(self):
        """
        Sets basic attributes directly dependent on "MONOMERS" dict from polymer
        config file.

        Raises:
            InvalidMonomerId: raised if monomer code does not match anything
                in polymer config file
        """

        #  init dict to store attrs and retrieve MONOMERS dict from polymer
        #  config
        basic_monomer_attrs = {}
        monomer_dict = self.__polymer_config["MONOMERS"]

        #  try loading monomer info from polymer config file, populate attrs
        #  dict
        try:
            basic_monomer_attrs["neutral_mass"] = monomer_dict[self.id][0]
            basic_monomer_attrs["func_groups"] = monomer_dict[self.id][1]
            basic_monomer_attrs["full_name"] = monomer_dict[self.id][2]
        except KeyError:
            raise InvalidMonomerId(
                monomer_id=self.id,
                available_monomers=monomer_dict.keys())

        #  set attrs in attrs dict
        for k, v in basic_monomer_attrs.items():
            setattr(self, k, v)

    def retrieve_loss_products(self):
        """
        Retrieve info on neutral losses for monomer

        Returns:
            List[float or str]: list of neutral losses associated with monomer,
                denoted either by functional string or mass.
        """
        if self.id in self.__polymer_config["LOSS_PRODUCTS"]:
            return self.__polymer_config[
                "LOSS_PRODUCTS"][self.id]
        return None

    def retrieve_sidechain_ions(self):
        """
        Check if monomer is ionizable at sidechain; if so return ionization info

        Returns:
            bool, List[str, float, float] or None
        """
        ionizable, ion_info = False, None
        if self.id in self.__polymer_config["IONIZABLE_SIDECHAINS"]:
            ion_info = self.__polymer_config[
                "IONIZABLE_SIDECHAINS"][self.id][self.__mode]
            if ion_info:
                ionizable = True

        return ionizable, ion_info

    def retrieve_ms2_signature_ions(self, signature_types):
        """
        Checks whether monomer has signature ion of specified type(s) and, if
        so, retrieves these signatures ions from polymer config file.

        Args:
            signature_types (List[str], optional): list of string identifiers
                for signature types

        Returns:
            List[float] (optional): list of associated signature
                ion m/z values for monomer, if applicable
        """
        if not signature_types:
            return None
        signatures = {}
        for sig_type in signature_types:
            if sig_type in self.__polymer_config[
                    "MS2_SIGNATURE_IONS"]:
                if self.id in self.__polymer_config[
                        "MS2_SIGNATURE_IONS"][sig_type]:
                    signatures[sig_type] = self.__polymer_config[
                        "MS2_SIGNATURE_IONS"][sig_type][self.id]
        return signatures

class Modification():
    def __init__(
        self,
        polymer_info,
        mod_id,
        mode
    ):
        """
        Defines covalent modification to be used in later silico functions.

        Args:
            polymer_info (dict): polymer-specific config dict.
            mod_id (str): string identifier for covalent modification
            mode (str): denotes mass spec mode - "pos" or "neg" for positive and
                negative, respectively.

        Props:
            id (str): modification id string (typically one letter code).
            mass (float): neutral monoisotopic mass for modification.
            termini (List[int]): list of termini modification can target.
            side_chain_attachments (List[str]): list of monomer id strings for
                monomers that modification can target @ sidechains.
            free_mod_fragments(List[float]): list of m/z values for dissociated
                modification MS2 fragments.
            mass_diff (float): mass lost upon addition of covalent modification
                defined as: (sequence_mass + modification_mass) -
                modified_sequence_mass.
            universal_ms2_shift (bool): specifies whether modification shifts
                every associated backbone ms2 fragment.
        """

        #  set id, mode and polymer_config props
        self.__mode = mode
        self.__polymer_config = load_polymer_info(polymer_info)
        self.id = mod_id

        #  retrieve relevant info on modification from polymer config, set
        #  config-dependent props
        mod_info = self.__polymer_config["MODIFICATIONS"][self.id]
        for k, v in mod_info.items():
            setattr(self, k, v)
