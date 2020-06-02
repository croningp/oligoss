"""
This file contains functions for the main Polymersoup experimental workflows
"""
import os
from ..silico.polymer_info.polymer import Polymer
from ..silico.ms1_silico import generate_ms1_ions
from ..silico.silico_handler import (
    get_ms2_silico_dict_from_compositions,
    combine_ms1_ms2_silico_dicts
)
from ..extractors.extractor_classes import RipperDict
from ..extractors.filters import (
    mzml_to_json,
    return_jsons,
    prefilter_all
)
from ..extractors.data_extraction import (
    extract_MS1_EICs,
    confirm_ms2_fragments_isomeric_sequences
)
from ..extractors.general_functions import write_to_json, open_json
from ..postprocessing.postprocess import postprocess_ripper

def exhaustive_screen(params, data_folder, out_folder):

    #  convert mzml files to ripper JSONs
    mzml_to_json(
        input_folder=data_folder,
        extractor_parameters=params.extractors
    )

    # #  get list of filepaths to ripper JSON files
    rippers = return_jsons(data_folder)

    filtered_rippers = [
        os.path.join(data_folder, "prefiltered", ripper_f)
        for ripper_f in os.listdir(os.path.join(data_folder, "prefiltered"))
    ]

    # pre-filter ripper JSONs and save output
    filtered_rippers = prefilter_all(
        rippers=rippers,
        extractor_parameters=params.extractors
    )

    #  get Polymer object from input parameters
    polymer = Polymer(params_obj=params)

    #  generate dict of compositions and MS1 precursors for screening
    compositional_ms1_dict = generate_ms1_ions(
        params=params,
        polymer=polymer,
        sequencing=False
    )

    ms1_silico_o_file = os.path.join(
        out_folder, "extracted", "MS1_pre_screening.json")
    if not os.path.exists(ms1_silico_o_file):
        os.mkdir(os.path.dirname(ms1_silico_o_file))

    #  write MS1 compositional silico dict to file
    write_to_json(
        write_dict=compositional_ms1_dict,
        output_json=ms1_silico_o_file
    )

    print(f"written MS1 dict to {ms1_silico_o_file}")
    #  iterate through each ripper file, screening for ms1 hits and then ms2
    #  fragments for potential hits
    for ripper in filtered_rippers:

        #  open ripper data, set output folder
        ripper_data = RipperDict(open_json(ripper))
        ripper_output = os.path.join(out_folder, os.path.basename(ripper))
        if not os.path.exists(ripper_output):
            os.mkdir(ripper_output)

        #  screem ripper for ms1 compositions
        ms1_hits = screen_ripper(
            ripper_data=ripper_data,
            compositional_ms1_dict=compositional_ms1_dict,
            polymer=polymer,
            params=params,
            ripper_output=ripper_output
        )

        #  screen ripper for ms2 fragments for compositions found at ms1
        screen_ripper_ms2(
            ms1_hits=ms1_hits,
            ripper_data=ripper_data,
            compositional_ms1_dict=compositional_ms1_dict,
            polymer=polymer,
            params=params,
            ripper_output_dir=ripper_output
        )

        postprocess_ripper(
            ripper_folder=ripper_output,
            postprocess_parameters=params.postprocess
        )

def screen_ripper(
    ripper_data,
    compositional_ms1_dict,
    polymer,
    params,
    ripper_output
):
    """
    Screens a single RipperDict object for MS1 precursors.

    Args:
        ripper_data (RipperDict): RipperDict object
        compositional_ms1_dict (Dict[str, List[float]]): dict of composition
            strings and corresponding precursor m/z values (MS1 ions)
        polymer (Polymer): Polymer object
        params (Parameters): Parameters object
        ripper_output (str): directory of output data folder for ripper

    Returns:
        List[str]: list of compositions that have been found in MS1 data of
            input ripper_data
    """
    #  get list of compositions with MS1 ions of sufficent abundance
    #  MS1 EICs for these sequences are also dumped to file
    ms1_hits = extract_MS1_EICs(
        MS1_silico=compositional_ms1_dict,
        extractor_parameters=params.extractors,
        MS1_dict=ripper_data.ms1,
        filename=os.path.basename(ripper_output),
        output_folder=ripper_output
    ).keys()

    return ms1_hits


def screen_ripper_ms2(
    ripper_data,
    ms1_hits,
    compositional_ms1_dict,
    ripper_output_dir,
    params,
    polymer
):
    """
    Screens single ripper object for MS1 precursor - MS2 product associations
    for sequences.

    Args:
        ripper_data (RipperDict): RipperDict object
        ms1_hits (List[str]): list of compositions that have been found in
            sufficient abundance in ripper_data MS1 spectra
        compositional_ms1_dict (Dict[str, List[float]]): dict of compositions
            and corresponding precursor m/z values
        ripper_output_dir (str): directory of output folder for ripper output
        params (Parameters): Parameters object
        polymer (Polymer): Polymer object
    """

    #  get silico dict of MS2 ions
    ms2_silico_dict = get_ms2_silico_dict_from_compositions(
        ms1_hits=ms1_hits,
        params=params,
        polymer=polymer
    )

    #  combine ms1 and ms2 silico dicts to create full MS/MS silico dict
    full_msms_silico_dict = combine_ms1_ms2_silico_dicts(
        ms1_silico_dict={
            k: v for k, v in compositional_ms1_dict.items()
            if k in ms1_hits
        },
        ms2_silico_dict=ms2_silico_dict
    )

    out_basename = os.path.basename(ripper_output_dir)

    #  set filename for MS/MS silico dict and dump to json
    out_silico_json = os.path.join(
        ripper_output_dir,
        f"{out_basename}_PRE_fragment_screening_silico_dict.json"
    )

    write_to_json(
        write_dict=full_msms_silico_dict,
        output_json=os.path.join(
            ripper_output_dir,
            out_silico_json
        )
    )

    #  init dict to store confirmed fragments and spectral matches for all
    #  sequences
    confirmed_fragment_dict, spectral_matches = {}, {}

    for composition in full_msms_silico_dict["compositions"]:

        confirmed_ms2 = confirm_ms2_fragments_isomeric_sequences(
            target_composition=composition,
            precursors=compositional_ms1_dict[composition],
            full_msms_silico_dict=full_msms_silico_dict,
            ripper_data=ripper_data,
            params=params
        )

        confirmed_fragment_dict.update(confirmed_ms2[0])
        spectral_matches.update(confirmed_ms2[1])

    write_to_json(
        write_dict=confirmed_fragment_dict,
        output_json=os.path.join(
            ripper_output_dir,
            f"{out_basename}_confirmed_fragment_dict.json"
        )
    )

    write_to_json(
        write_dict=spectral_matches,
        output_json=os.path.join(
            ripper_output_dir,
            f"{out_basename}_spectra_matches.json"
        )
    )
