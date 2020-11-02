"""
This file contains tests for Polymer class

"""
import os
import re
import copy
import pytest

from oligoss.utils.parameter_handlers.parameter_handlers import generate_parameters
from oligoss.utils.parameter_handlers.polymer_param_handlers import load_polymer_info
from oligoss.utils.errors import InvalidMonomerId, InvalidModificationTarget
from oligoss.silico.polymer_info.polymer import Polymer, Monomer, Modification

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
def test_invalid_polymer_alias(simple_params):
    """
    Tests for invalid polymer alias being picked up and Exception raised.
    """

    #  set up params with undefined polymer class and ensure appropriate error
    #  is raised
    fake_params = copy.deepcopy(simple_params)
    fake_params.polymer_class = "DNA"

    with pytest.raises(Exception, match=r"valid polymer"):
        err_polymer = Polymer(fake_params)
        for key in err_polymer.__dict__:
            assert key == "params"

@pytest.mark.unit
def test_monomers_match_parameters(simple_params):
    """
    Test to make sure monomers match between input parameters and Polymers
    object

    Args:
        simple_params (Parameters): Parameters object
    """
    polymer = Polymer(simple_params)

    assert not polymer.modifications

    monomer_ids = [monomer for monomer in simple_params.monomers]
    polymer_monomer_ids = [
        monomer.id for monomer in polymer.monomers]
    assert sorted(monomer_ids) == sorted(polymer_monomer_ids)

    for monomer in polymer.monomers:
        copy_monomer = Monomer(
            monomer_id=monomer.id,
            mode=simple_params.mode,
            ms2_signature_types=simple_params.silico.ms2.signatures,
            polymer_info=simple_params.polymer_class,
            modification_targets=None)

        for prop, value in monomer.__dict__.items():
            if prop.find("__") == -1:
                if copy_monomer.__dict__[prop] != value:
                    raise Exception(prop)


@pytest.mark.unit
def test_invalid_monomers(simple_params):
    """
    Test to ensure appropriate error is raised if invalid monomer ids are
    supplied in parameters.

    Args:
        simple_params (Parameters): Parameters object
    """

    #  make copy of simple_params for changing monomers
    params = copy.deepcopy(simple_params)

    #  random list of invalid monomers not pre-configured
    invalid_monomers = [
        ["B", "X", "Y"],
        ["A", "D", "UEURUA"],
        ["2", "62", 21],
        ["4", "4", "xylophone"]
    ]

    #  assert that each list of invalid monomers raised InvalidMonomerId custom
    #  error
    for invalid_list in invalid_monomers:
        with pytest.raises(
            InvalidMonomerId,
            match=r"Here are the valid monomer one letter codes"
        ):
            params.monomers = invalid_list
            Polymer(params_obj=params)

@pytest.mark.unit
def test_invalid_modification_target(simple_params, polymer_info):
    """
    Test to ensure invalid modification targets are handled appropriately - if
    a proposed modification and target do not match, raise
    InvalidModificationTarget error.

    Args:
        simple_params (Parameters): Parameters object
        polymer_info (dict): depsipeptide polymer config dict
    """
    #  make copy of simple params for changing modifications
    params = copy.deepcopy(simple_params)

    #  list of invalid sidechain modification targets
    invalid_sidechain_targets = [
        {
            "E": ["Pal", "Ole", "Ace"],
            "D": ["Pal", "Ole", "Ace"]
        },
        {
            "R": ["Fmc", "Bnz", "BtA"],
            "K": ["Fmc", "Bnz", "BtA"]
        },
        {
            "A": ["Pal", "Ole", "Ace", "Fmc", "Bnz", "BtA"]
        }
    ]

    #  list of invalid terminal modification targets
    invalid_terminal_targets = [
        {
            0: ["Fmc", "BtA"]
        },
        {
            -1: ["Pal", "Ole", "Ace"]
        },
        {
            0: ["Fmc", "BtA"],
            -1: ["Pal", "Ole", "Ace"]
        }
    ]

    #  list of valid modification targets
    valid_targets = [
        {
            "E": ["Fmc", "Bnz", "BtA"],
            "D": ["Fmc", "Bnz", "BtA"]
        },
        {
            "K": ["Pal", "Ole", "Ace"],
            "R": ["Pal", "Ole", "Ace"]
        },
        {
            0: ["Pal", "Ole", "Ace"],
            -1: ["Fmc", "Bnz", "BtA"]
        },
        {
            "E": ["Fmc", "Bnz", "BtA"],
            "D": ["Fmc", "Bnz", "BtA"],
            "K": ["Pal", "Ole", "Ace"],
            "R": ["Pal", "Ole", "Ace"],
            0: ["Pal", "Ole", "Ace"],
            -1: ["Fmc", "Bnz", "BtA"]
        }
    ]

    #  combine targets
    all_targets = (invalid_terminal_targets + invalid_sidechain_targets
                   + valid_targets)

    #  iterate through target dicts and make sure each are handled properly
    for target_dict in all_targets:

        #  set modification parameters to current targets
        params.silico.modifications = target_dict

        #  ensure Polymer instantiated without issue if targets are valid
        if target_dict in valid_targets:
            polymer = Polymer(params_obj=params)

            #  check that polymer modifications have been instantiated
            assert polymer.modifications

            #  iterate through valid modifications and check that they target
            #  the correct places
            for target, modifications in polymer.modifications.items():

                #  check that each modification is a Modification object and
                #  is targeted appropriately
                for modification in modifications:
                    assert type(modification) == Modification

                    modification.termini = [
                        str(x) for x in modification.termini]

                    if target not in modification.termini:
                        assert target in modification.side_chain_attachments
                    else:
                        assert target not in modification.side_chain_attachments

                    if target not in modification.side_chain_attachments:
                        assert target in modification.termini
                    else:
                        assert target not in modification.termini

                    #  make sure Modification attributes match what has been
                    #  specified in polymer-specific config file
                    for prop, value in modification.__dict__.items():
                        if prop != "id" and not re.search("Modification", prop):
                            if prop == "termini":
                                value = [int(x) for x in value]
                            assert value == polymer_info[
                                "MODIFICATIONS"][modification.id][prop]

        #  ensure appropriate error is raised if invalid terminal targets are
        #  provided
        if target_dict in invalid_terminal_targets:

            with pytest.raises(
                    InvalidModificationTarget):
                polymer = Polymer(params_obj=params)

        #  ensure appropriate error is raised if invalid sidechain targets are
        #  provided
        if target_dict in invalid_sidechain_targets:

            with pytest.raises(
                    InvalidModificationTarget,
                    match=r"sidechain"):
                polymer = Polymer(params_obj=params)
