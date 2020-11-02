"""
This file contains tests for Monomer class
"""

import os
import pytest

from oligoss.utils.parameter_handlers.parameter_handlers import generate_parameters
from oligoss.utils.parameter_handlers.polymer_param_handlers import load_polymer_info
from oligoss.silico.polymer_info.polymer import Monomer

HERE = os.path.abspath(os.path.dirname(__file__))
INPUTS_FOLDER = os.path.join(HERE, '..', 'input_param_files')

@pytest.fixture
def simple_params():
    return generate_parameters(
        params_json=os.path.join(
            INPUTS_FOLDER,
            'full_input_parameters_no_defaults.json'))

@pytest.fixture
def polymer_info():
    return load_polymer_info(polymer_class="depsipeptide")

@pytest.mark.unit
def test_monomer_props(simple_params, polymer_info):
    """
    Test for simple monomer props

    Args:
        simple_params (Parameters): Parameters object
        polymer_info (dict): polymer info dict from config file
    """
    for monomer_id in simple_params.monomers:
        monomer = Monomer(
            monomer_id=monomer_id,
            polymer_info=simple_params.polymer_class,
            mode=simple_params.mode,
            ms2_signature_types=simple_params.silico.ms2.signatures,
            modification_targets=None
        )
        assert monomer.id == monomer_id
        assert monomer.neutral_mass == polymer_info["MONOMERS"][monomer.id][0]
        assert monomer.func_groups == polymer_info["MONOMERS"][monomer.id][1]
        assert monomer.full_name == polymer_info["MONOMERS"][monomer.id][2]

@pytest.mark.unit
def test_monomer_neutral_loss(simple_params, polymer_info):
    """
    Test for monomer neutral loss products

    Args:
        simple_params (Parameters): Parameters object
        polymer_info (dict): polymer info dict from config file

    """
    #  list of monomers with no neutral losses
    test_monomers = [
        ["A", "V", "G"],
        ["L", "I", "F"],
        ["W", "Y", "L", "V", "G", "g"]
    ]

    #  list of monomers with neutral losses
    loss_monomers = [
        ["S", "T", "E", "N", "Q"],
        ["K", "S", "R"]
    ]

    #  add loss monomers to test monomers
    test_monomers += loss_monomers

    #  cheeky wee internal function to generate Monomer objects from id
    def generate_monomer_obj_from_id(id):
        return Monomer(
            monomer_id=id,
            polymer_info=polymer_info,
            mode=simple_params.mode,
            ms2_signature_types=simple_params.silico.ms2.signatures,
            modification_targets=None)

    #  iterate through monomer tests and make sure loss products are correct
    for monomer_list in test_monomers:
        monomers = map(
            lambda x: generate_monomer_obj_from_id(x), monomer_list)
        for monomer in monomers:
            if monomer_list not in loss_monomers:
                assert not monomer.loss_products
                assert monomer.id not in polymer_info[
                    "LOSS_PRODUCTS"]
            else:
                if not monomer.loss_products:
                    raise Exception(monomer.id)
                assert monomer.loss_products
                assert monomer.loss_products == polymer_info[
                    "LOSS_PRODUCTS"][monomer.id]

@pytest.mark.unit
def test_ionizable_sidechains(simple_params, polymer_info):
    """
    Test for ionizable sidechain info in positive, negative mode.

    Args:
        simple_params (Parameters): Parameters object
        polymer_info (dict): polymer info dict from config file
    """
    #  list of monomers that do not have ionizable sidechains
    neutral_monomers = [
        "A", "V", "G", "L", "I"
    ]

    #  list of monomers that ionize in positive mode
    cationic_monomers = [
        "K", "R", "H"
    ]

    #  list of monomers that ionize in negative mode
    anionic_monomers = [
        "D", "E"
    ]

    #  combine monomer lists
    all_monomers = neutral_monomers + cationic_monomers + anionic_monomers

    #  check that ionization for monomers is correct in both positive and
    #  negative mode
    for mode in ["pos", "neg"]:
        for monomer in all_monomers:
            monomer_obj = Monomer(
                monomer_id=monomer,
                polymer_info=polymer_info,
                mode=mode,
                ms2_signature_types=simple_params.silico.ms2.signatures,
                modification_targets=None)
        if monomer in neutral_monomers:
            assert not monomer_obj.ionizable
            assert not monomer_obj.ion_info
        elif monomer in cationic_monomers:
            if mode == "neg":
                assert not monomer_obj.ionizable
                assert not monomer_obj.ion_info
            else:
                assert monomer_obj.ionizable
                assert monomer_obj.ion_info == polymer_info[
                    "IONIZABLE_SIDECHAINS"][monomer_obj.id][mode]
        elif monomer in anionic_monomers:
            if mode == "pos":
                assert not monomer_obj.ionizable
                assert not monomer_obj.ion_info
            else:
                assert monomer_obj.ionizable
                assert monomer_obj.ion_info == polymer_info[
                    "IONIZABLE_SIDECHAINS"][monomer_obj.id][mode]

@pytest.mark.unit
def test_invalid_monomers(simple_params, polymer_info):
    pass
