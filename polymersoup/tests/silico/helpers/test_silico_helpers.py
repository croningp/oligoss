"""
This file contains tests for silico module helper functions. These are basic
functions required to perform simple tasks in the silico module, such as
calculate sequence masses or manipulate sequence strings.
"""

import os
import copy
import pytest
import random

from ....utils.parameter_handlers.parameter_handlers import generate_parameters
from ....utils.parameter_handlers.polymer_param_handlers import (
    load_polymer_info
)
from ....utils.global_chemical_constants import FUNCTIONAL_GROUPS, CATIONS
from ....silico.polymer_info.polymer import Polymer
from ....silico.helpers.helpers import (
    find_sequence_mass,
    return_core_sequence,
    reverse_sequence,
    ionize_sequence_precursors,
    add_sidechain_neutral_loss_products_sequence,
    generate_all_sequences
)

HERE = os.path.abspath(os.path.dirname(__file__))
INPUTS_FOLDER = os.path.join(os.path.dirname(HERE), '..', 'input_param_files')

@pytest.fixture
def params():
    return generate_parameters(
        params_json=os.path.join(
            INPUTS_FOLDER,
            'full_input_parameters_with_silico_modifications.json'))

@pytest.fixture
def test_sequences():
    """
    Dict of depsipeptide sequence strings and associated neutral monoisotopic
    masses for testing silico helper functions. Masses were calculated from:
    http://db.systemsbiology.net:8080/proteomicsToolkit/FragIonServlet.html

    Returns:
        Dict[str, float]: dict of sequence strings and masses
    """
    return {
        "GLESGA": 532.24934,
        "MARYHILL": 1015.56375,
        "WESTEND": 879.32468,
        "LANARKSHIRE": 1293.72661,
        "MANCHESTER": 1176.46524,
        "GHANA": 468.20814,
        "NIGERIA": 771.42394,
        "YEMEN": 684.24254,
        "ALGERIA": 728.41813,
        "EGYPT": 565.23844,
        "ENGLAND": 731.30864,
        "HAITI": 553.32244,
        "IRAQ": 486.29147,
        "IRAN": 472.27582,
        "WALES": 604.28572,
        "IRELAND": 829.42942,
        "AMERICA": 792.36227,
        "CANADA": 563.20101,
        "GERMANY": 839.35963,
        "INDIA": 544.28572,
        "SPAIN": 500.25951,
        "AFGHANISTAN": 1101.52036
    }

@pytest.fixture
def modified_test_sequences():
    """
    Modified sequence strings for testing.

    Returns:
        Dict[str, List[str]]: keys = modified sequence strings. value[0] =
        unmodified core sequence, value[1] = reversed modified sequence string
    """
    return {
        "[Pal]GLE(Fmc)S(Trt)GA[Fmc]": [
            "GLESGA",
            "[Fmc]AGS(Trt)E(Fmc)LG[Pal]"
        ],
        "[Ole]MARYHILL[Bnz]": [
            "MARYHILL",
            "[Bnz]LLIHYRAM[Ole]"
        ],
        "[Ace]WE(Bnz)S(Trt)T(Trt)E(Fmc)ND": [
            "WESTEND",
            "DNE(Fmc)T(Trt)S(Trt)E(Bnz)W[Ace]"
        ],
        "[Boc]LANAR(Tfa)KSHIR(Tfa)E": [
            "LANARKSHIRE",
            "ER(Tfa)IHSKR(Tfa)ANAL[Boc]"
        ],
        "[Ole]MANCHE(Bnz)S(Trt)T(Trt)ER(Tfa)[Fmc]": [
            "MANCHESTER",
            "[Fmc]R(Tfa)ET(Trt)S(Trt)E(Bnz)HCNAM[Ole]"
        ],
        "CIR(Tfa)E(BtA)NCE(Bnz)S(Trt)T(Trt)E(BtA)R(Ace)[BtA]": [
            "CIRENCESTER",
            "[BtA]R(Ace)E(BtA)T(Trt)S(Trt)E(Bnz)CNE(BtA)R(Tfa)IC"
        ]
    }

@pytest.fixture
def polymer(params, test_sequences):
    params.monomers = polymer.monomers = list(
        set("".join(test_sequences.keys())))
    test_polymer = Polymer(params_obj=params)
    if type(test_polymer.mass_diff) == str:
        test_polymer.mass_diff = FUNCTIONAL_GROUPS[test_polymer.mass_diff]
    return test_polymer

@pytest.mark.unit
def test_sequence_masses(test_sequences, polymer):
    """
    Tests find sequence mass function against pre-defined peptide sequences with
    known neutral monoisotopic masses.

    Args:
        test_sequences (Dict[str, float]): test_sequences fixture
        polymer (Polymer): polymer fixture
    """

    #  iterate through test sequences and check correct neutral monoisotopic
    #  mass is calculated to 3 decimal places - NOTE: this needs to be improved
    #  to 4 or 5 decimal places ASAP
    for sequence in test_sequences:

        sequence_mass = find_sequence_mass(
            sequence=sequence,
            polymer=polymer,
            n_rounded=3)

        assert sequence_mass == round(test_sequences[sequence], 3)

@pytest.mark.unit
def test_core_sequence(modified_test_sequences):
    """
    Tests return_core_sequence function, which should remove all terminal
    and sidechain modifications from a sequence string and return the "core"
    sequence comprised solely of the monomer one-letter codes.

    Args:
        modified_test_sequences (Dict[str, List[str]]): modified test sequences
            fixture
    """
    for modified_seq, unmodified_seq in modified_test_sequences.items():
        assert return_core_sequence(modified_seq) == unmodified_seq[0]

@pytest.mark.unit
def test_reverse_sequence(modified_test_sequences):
    """
    Tests reverse sequence function for modified sequences - codes for sidechain
    and terminal modifications must be correct in reversed sequence

    Args:
        modified_test_sequences (Dict[str, List[str]]): modified_test_sequences
            fixture.
    """
    for sequence, targets in modified_test_sequences.items():
        assert reverse_sequence(sequence) == targets[1]

@pytest.mark.unit
def test_ionize_sequence(test_sequences, polymer, params):
    """
    Tests ionize sequence function. This function should produce a full list
    of MS1 m/z values for target sequence, including adducts and neutral losses.

    Args:
        test_sequences (Dict[str, float]): test_sequences fixture
        polymer (Polymer): Polymer object
        params (Parameters): Parameters object
    """
    for sequence in test_sequences:
        mz_list = ionize_sequence_precursors(
            sequence=sequence,
            params=params,
            polymer=polymer
        )
        if test_sequences[sequence] < 1000:
            assert round(mz_list[0], 3) == round(
                test_sequences[sequence] + CATIONS["Na"][0], 3)
        else:
            assert round(mz_list[0], 2) == round(
                test_sequences[sequence] + CATIONS["Na"][0], 2)

@pytest.mark.unit
def test_neutral_losses(params):
    """
    Tests add_sidechain_neutral_loss_products_sequence function. This should
    take a sequence string, corresponding full neutral mass, and work out
    a list of neutral loss products from constituent monomers (including
    sidechain mods) and Parameters object.

    Args:
        params (Parameters): Parameters object
    """
    #  make copy of params for editing
    test_params = copy.deepcopy(params)

    #  increase neutral loss limit
    test_params.silico.ms1.max_neutral_losses = 6

    #  set modifications to target test sequences with sidechain mods
    test_params.silico.modifications = {
        "K": ["Ole", "Pal", "Ace"],
        "R": ["Ole", "Pal", "Ace"],
        "S": ["Trt"],
        "T": ["Trt"],
        "E": ["Bnz", "BtA"],
        "D": ["Bnz", "BtA"]
    }

    #  keys = sequences prone to neutral losses
    neutral_loss_seqs = {
        "KARKSN": [
            5,  # expected number of neutral loss products for "KARKSN"
            {
                "K(Ole)AR(Pal)S(Trt)N": 1,  # modified sequences, expected
                "K(Pal)AR(Ole)S(Trt)N": 1,  # number of neutral losses
                "K(Ace)AR(Pal)S(Trt)N": 1,
                "K(Ole)AR(Pal)SN": 2,
                "KAR(Ole)S(Trt)N": 2,
                "K(Ace)ARS(Trt)N": 2,
                "K(Ole)ARSN": 3,
                "KARS(Trt)N": 3,
                "K(Ace)ARSN": 3,
                "KARSN": 4
            }
        ],
        "RTEDSN": [
            6,
            {
                "R(Ole)T(Trt)E(Bnz)D(BtA)S(Trt)N": 1,
                "R(Ace)T(Trt)E(BtA)D(Bnz)S(Trt)N": 1,
                "R(Pal)T(Trt)E(Bnz)D(BtA)S(Trt)N": 1,
                "R(Ole)T(Trt)E(Bnz)D(BtA)SN": 2,
                "RT(Trt)E(BtA)D(Bnz)S(Trt)N": 2,
                "R(Pal)T(Trt)ED(BtA)S(Trt)N": 2,
                "RT(Trt)E(Bnz)D(BtA)SN": 3,
                "RT(Trt)E(BtA)D(Bnz)SN": 3,
                "R(Pal)T(Trt)EDS(Trt)N": 3,
                "RT(Trt)EDS(Trt)N": 4,
                "RTE(BtA)D(Bnz)SN": 4,
                "R(Pal)T(Trt)EDSN": 4,
                "RTEDS(Trt)N": 5,
                "RTE(BtA)DSN": 5,
                "R(Pal)TEDSN": 5,
                "RTEDSN": 6,
            }
        ]
    }

    #  make sure monomers list contains every monomer used in tests
    test_params.monomers = list(set(
        "".join(test_params.silico.modifications.keys())
        + "".join(neutral_loss_seqs.keys())))

    #  Polymer object using new params
    test_polymer = Polymer(params_obj=test_params)

    #  iterate through modified and unmodified sequences. Check that all
    #  neutral losses are present for unmodified seqs. Modified seqs should
    #  only have two neutral losses - OH, H2O for asparagine (N)
    for seq, info in neutral_loss_seqs.items():
        full_neutral_mass = find_sequence_mass(
            sequence=seq,
            polymer=test_polymer)
        neutral_losses = add_sidechain_neutral_loss_products_sequence(
            sequence=seq,
            sequence_mass=full_neutral_mass,
            polymer=test_polymer,
            params=test_params
        )

        assert len(neutral_losses) == info[0]

        for mod_seq in info[1]:
            mod_neutral_mass = find_sequence_mass(
                sequence=mod_seq,
                polymer=test_polymer
            )
            mod_neutral_masses = add_sidechain_neutral_loss_products_sequence(
                sequence=mod_seq,
                sequence_mass=mod_neutral_mass,
                polymer=test_polymer,
                params=test_params
            )
            if len(mod_neutral_masses) != info[1][mod_seq]:
                raise Exception(mod_seq, len(
                    mod_neutral_masses), info[1][mod_seq])
            assert len(mod_neutral_masses) == info[1][mod_seq]

@pytest.mark.unit
def test_generate_all_sequence(params):
    """
    Tests generate_all_sequences function

    Args:
        params (Parameters): Parameters object
    """
    #  make copy of params for editing
    local_params = copy.deepcopy(params)
    local_params.silico.min_length = 1

    #  load polymer info and choose 10 monomers at random
    polymer_info = load_polymer_info(
        polymer_class=params.polymer_class)

    #  set 10 random monomers in params and polymer
    local_params.monomers = random.sample(polymer_info[
        "MONOMERS"].keys(), k=7)
    local_params.silico.isomeric_targets = None
    local_polymer = Polymer(params_obj=local_params)

    def check_length_distribution():
        N_sequence_space = 0
        for x in range(1, local_params.silico.max_length + 1):
            N_sequence_space += (len(local_params.monomers)
                                 ** x)

        return N_sequence_space

    for x in range(1, len(local_params.monomers)):
        local_params.silico.max_length = x
        sequences = generate_all_sequences(
            polymer=local_polymer,
            params=local_params,
            sequencing=True
        )

        assert len(sequences) == check_length_distribution()
