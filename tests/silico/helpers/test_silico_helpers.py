"""
This file contains tests for silico module helper functions. These are basic
functions required to perform simple tasks in the silico module, such as
calculate sequence masses or manipulate sequence strings.
"""

import os
import copy
import pytest
import random

from oligoss.utils.parameter_handlers.parameter_handlers import generate_parameters
from oligoss.utils.parameter_handlers.polymer_param_handlers import (
    load_polymer_info
)
from oligoss.utils.global_chemical_constants import FUNCTIONAL_GROUPS, CATIONS
from oligoss.silico.polymer_info.polymer import Polymer
from oligoss.silico.helpers.helpers import (
    find_sequence_mass,
    return_core_sequence,
    reverse_sequence,
    ionize_sequence_precursors,
    generate_all_sequences,
    get_isomeric_seqs,
    get_full_neutral_masses_sequence,
    add_terminal_mods_strings,
    get_composition,
    return_monomer_ids_sequence
)

HERE = os.path.abspath(os.path.dirname(__file__))
INPUTS_FOLDER = os.path.join(os.path.dirname(HERE), '..', 'input_param_files')

@pytest.fixture
def params():
    return generate_parameters(
        params_json=os.path.join(
            INPUTS_FOLDER,
            'full_input_parameters_without_silico_modifications.json'))

@pytest.fixture
def mods_params():
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
def test_sequence_list_with_sidechain_mods(mods_params):
    """
    Tests generation of sequence strings from monomers with targeted sidechain
    modifications.

    Args:
        mods_params (Parameters): parameters object with targeted sidechain
            modifications.
    """
    #  instantiate Polymer object with modifications in Parameters object
    mod_polymer = Polymer(params_obj=mods_params)

    #  generate full list of sequence strings with modified sidechains included
    mod_seq_pool = generate_all_sequences(
        polymer=mod_polymer,
        params=mods_params,
        sequencing=True)

    #  get monomer ids as stored in Polymer object
    monomer_ids = [
        monomer.id for monomer in mod_polymer.monomers
    ]

    #  check monomer id strings contain appropriate modified monomer ids
    assert sorted(monomer_ids) == sorted(
        ["R", "N", "K(Ole)", "K(Pal)", "K(Ace)", "S(Trt)", "E(Bnz)"])

    #  iterate through sequence strings and get unique monmer id strings
    for seq in mod_seq_pool:
        monomers = return_monomer_ids_sequence(
            sequence=seq,
            return_modified=True,
            return_set=True
        )

        #  check that monomers targeted for modification are not present in
        #  unmodified form
        assert not [
            monomer for monomer in monomers
            if monomer in ["K", "S", "E"]
        ]

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
def test_full_neutral_masses(params, polymer):
    """
    Test for neutral masses. Checks whether function
    "get_full_neutral_masses_sequence" calculates accurate full monoisotopic
    masses and sidechain-specific loss products. Also tests whether
    max_neutral_losses is enforced from silico parameters.

    Args:
        params (Parameters): Parameters object
        polymer (Polymer): Polymer object
    """
    #  make a local copy of Parameters object for changing parameters
    local_params = copy.deepcopy(params)
    local_params.monomers = ["A", "S", "G", "N", "Q", "V"]

    #  set upper limit on neutral loss cap to apply in local_params
    max_loss_caps = 10

    #  dict of sequence strings and associated full neutral masses, including
    #  all possible loss products
    loss_prone_sequences = {
        "AAV": [259.15322],
        "ASGNQ": [
            475.20273,
            458.17618,
            457.192165,
            441.14963,
            440.165615,
            439.1816,
            423.139065,
            422.15505
        ]
    }

    #  iterate through sequenes and associated neutral masses, ensure number
    #  of neutral sequences = full neutral mass + all permissible loss products
    for seq, neutral_masses in loss_prone_sequences.items():
        for i in range(max_loss_caps):
            local_params.silico.ms1.max_neutral_losses = i
            test_masses = get_full_neutral_masses_sequence(
                sequence=seq,
                polymer=Polymer(local_params),
                params=local_params,
                ms_level=1
            )

            #  losses should always be capped by parameters override
            assert len(test_masses) == min(i + 1, len(neutral_masses))

            #  ensure all neutral masses are accounted for - no nasty surprises
            assert not [
                mass for mass in test_masses
                if mass not in neutral_masses]

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
            N_sequence_space += (len(local_params.monomers) ** x)

        return N_sequence_space

    for x in range(1, len(local_params.monomers)):
        local_params.silico.max_length = x
        sequences = generate_all_sequences(
            polymer=local_polymer,
            params=local_params
        )
        sequences = sorted([seq for seq in sequences])
        assert sequences == sorted(list(set(sequences)))
        assert len(sequences) == check_length_distribution()

@pytest.mark.unit
def test_isomeric_targets(params):
    """
    Test for isomeric sequence lists.

    Args:
        params (Parameters): Parameters object.
    """

    local_params = copy.deepcopy(params)

    #  list of target sequence lists
    targets = [
        ["AAAA", "TRED"],
        ["NVTEDPQ", "GAGA", "GAGA", "GAGA"],
        ["KICKERHV"]
    ]

    #  iterate through isomeric targets and check get_isomeric_seqs for each
    for test_data in targets:

        local_params.silico.isomeric_targets = test_data

        #  generate list of sequences and compositions that are isomeric to
        #  one or more targets
        isomeric_seqs = get_isomeric_seqs(
            target_sequences=test_data,
            polymer=polymer
        )
        isomeric_compositions = [comp for comp in generate_all_sequences(
            polymer=polymer,
            params=local_params,
            sequencing=False
        )]

        #  get composition strings that match isomers in test data
        compositions = list(set([
            "".join(sorted(target_seq)) for target_seq in test_data
        ]))

        #  check all isomeric seqs are  correct
        for seq in isomeric_seqs:
            assert "".join(sorted(seq)) in compositions

        #  check function works for composition strings
        assert sorted(compositions) == sorted(isomeric_compositions)

@pytest.mark.unit
def test_add_terminal_mods_strings():

    targets = ["FLOP", "F", "PLOP"]

    terminal_modifications = {"0": ["Yas", "Wow"], "-1": "Boc"}

    universal_positive_test = sorted(add_terminal_mods_strings(
        sequence_list=targets,
        terminal_modifications=terminal_modifications,
        universal_mod=True))

    universal_positive = sorted(
        ["[Wow]FLOP[Boc]", "[Wow]F[Boc]", "[Wow]PLOP[Boc]", "[Yas]FLOP[Boc]",
            "[Yas]F[Boc]", "[Yas]PLOP[Boc]"])

    non_universal_positive_test = sorted(add_terminal_mods_strings(
        sequence_list=targets,
        terminal_modifications=terminal_modifications,
        universal_mod=False))

    non_universal_positive = sorted(
        ["[Yas]FLOP[Boc]", "[Wow]FLOP[Boc]", "[Wow]F[Boc]", "[Yas]PLOP",
            "[Wow]F", "F[Boc]", "FLOP", "[Yas]FLOP", "F", "[Yas]PLOP[Boc]",
            "[Yas]F", "PLOP", "[Yas]F[Boc]", "[Wow]PLOP[Boc]", "[Wow]PLOP",
            "PLOP[Boc]", "[Wow]FLOP", "FLOP[Boc]"])

    negative_test = sorted(add_terminal_mods_strings(
        sequence_list=targets,
        terminal_modifications={},
        universal_mod=True))

    assert universal_positive_test == universal_positive
    assert non_universal_positive_test == non_universal_positive
    assert negative_test == sorted(targets)

@pytest.mark.unit
def test_get_composition():

    # monomers should appear alphabetically
    test_sidechain = get_composition(sequence='AB(Pop)VER')

    # terminally modified monomers will remain in that position
    # other monomers should be in alphabetical order
    test_terminal = get_composition(sequence='[Boc]UTENEC[Yes]')
    test_double_mod = get_composition(sequence='AB(Pop)VER[Trt]')
    negative_test = get_composition(sequence='BDAEC')

    assert test_sidechain == 'AB(Pop)ERV'
    assert test_terminal == '[Boc]UEENTC[Yes]'
    assert test_double_mod == 'AB(Pop)EVR[Trt]'
    assert negative_test == 'ABCDE'
