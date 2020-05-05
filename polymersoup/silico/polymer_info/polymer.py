"""
This file contains classes for storing information on polymer-specific
features that are passed in from parameter handlers and passed on to
silico modules
"""
from ...utils.parameter_handlers.polymer_param_handlers import load_polymer_info
from ...utils.errors import InvalidMonomerId

class Polymer():
    def __init__(
        self,
        params_obj
    ):
        self.params = params_obj
        self.polymer_config = load_polymer_info(params_obj.polymer_class)
        self.monomers = self.instantiate_monomers()
        self.dominant_signatures = self.retrieve_dominant_signatures()
        self.mass_diff = self.polymer_config["MASS_DIFF"]
        self.elongation_unit = self.polymer_config["ELONGATION"]
        self.symmetry = self.polymer_config["SYMMETRY"]
        self.chain_terminators = self.polymer_config["CHAIN_TERMINATORS"]
        self.fragment_info = self.polymer_config["FRAG_SERIES"]

        self.sidechain_modifications, self.terminal_modifications = (
            self.instantiate_modifications())

    def instantiate_monomers(self):
        """
        Returns monomers parameter for instance of Polymer class

        Returns:
            List[Monomer]: list of Monomer objects
        """
        monomer_ids = self.params.monomers
        monomers = [
            Monomer(
                monomer_id=monomer_id,
                polymer_info=self.polymer_config,
                mode=self.params.mode,
                ms2_signature_types=self.params.silico.ms2.signatures)
            for monomer_id in monomer_ids]
        return monomers

    def retrieve_dominant_signatures(self):
        """
        Get list of monomer str codes for monomers expected to have dominant
        ms2 signatures

        Returns:
            List[str]: list of monomer str codes
        """
        dominant_sigs = []
        if not self.params.silico.ms2.signatures:
            return None

        for sig in self.params.silico.ms2.signatures:
            for monomer in self.monomers:
                if monomer.id in self.polymer_config[
                        "MS_SIGNATURE_IONS"][sig]["dominant"]:
                    dominant_sigs.append(monomer.id)
        return dominant_sigs

    def instantiate_modifications(self):
        """
        Reads params and passes on list of instances of Modifcation class for
        sidechain modifications and terminal modifications

        Returns:
            List[Modification], List[Modification]: sidechain_modifications,
                terminal_modifications
        """
        sidechain_mods, terminal_mods = [], []

        sidechain_targets = self.params.silico.ms1.sidechain_modifications
        terminal_targets = self.params.silico.ms1.terminal_modifications

        if sidechain_targets:
            sidechain_mods.extend(
                [
                    Modification(
                        polymer_info=self.polymer_config,
                        mod_id=mod,
                        mode=self.params.mode
                    )
                    for mod in sidechain_targets])
        if terminal_targets:
            terminal_mods.extend([
                Modification(
                    polymer_info=self.polymer_config,
                    mod_id=mod,
                    mode=self.params.mode
                )
                for mod in terminal_targets])

        return sidechain_mods, terminal_mods

class Monomer():
    def __init__(
        self,
        monomer_id,
        polymer_info,
        mode,
        ms2_signature_types
    ):
        if type(polymer_info) == str:
            self.polymer_config = load_polymer_info(polymer_info)
        else:
            self.polymer_config = polymer_info
        self.check_polymer_config
        self.mode = mode
        self.id = monomer_id
        self.set_basic_attrs()
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
        basic_monomer_attrs = {}
        monomer_dict = self.polymer_config["MONOMERS"]

        try:
            basic_monomer_attrs["neutral_mass"] = monomer_dict[self.id][0]
            basic_monomer_attrs["func_groups"] = monomer_dict[self.id][1]
            basic_monomer_attrs["full_name"] = monomer_dict[self.id][2]
        except KeyError:
            raise InvalidMonomerId(
                monomer_id=self.id,
                available_monomers=monomer_dict.keys())
        for k, v in basic_monomer_attrs.items():
            setattr(self, k, v)

    def check_polymer_config(self):
        if type(self.polymer_config) == str:
            self.polymer_config = load_polymer_info(self.polymer_config)

    def retrieve_loss_products(self):

        if self.id in self.polymer_config["LOSS_PRODUCTS"]:
            return self.polymer_config[
                "LOSS_PRODUCTS"][self.id]
        return None

    def retrieve_sidechain_ions(self):

        ionizable, ion_info = False, None
        if self.id in self.polymer_config["IONIZABLE_SIDECHAINS"]:
            ion_info = self.polymer_config[
                "IONIZABLE_SIDECHAINS"][self.id][self.mode]
            if ion_info:
                ionizable = True

        return ionizable, ion_info

    def retrieve_ms2_signature_ions(self, signature_types):
        if not signature_types:
            return None
        signatures = {}
        for sig_type in signature_types:
            if sig_type in self.polymer_config[
                    "MS2_SIGNATURE_IONS"]:
                if self.id in self.polymer_config[
                        "MS2_SIGNATURE_IONS"][sig_type]:
                    signatures[sig_type] = self.polymer_config[
                        "MS2_SIGNATURE_IONS"][sig_type][self.id]
        return signatures

class Modification():
    def __init__(
        self,
        polymer_info,
        mod_id,
        mode
    ):

        self.mode = mode
        self.polymer_config = load_polymer_info(polymer_info)
        self.id = mod_id

        mod_info = self.polymer_config["MODIFICATIONS"][self.id]
        for k, v in mod_info.items():
            setattr(self, k, v)
