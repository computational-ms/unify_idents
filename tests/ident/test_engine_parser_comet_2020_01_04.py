#!/usr/bin/env python
from pathlib import Path

from unify_idents.engine_parsers.ident.comet_2020_01_4_parser import (
    Comet_2020_01_4_Parser,
)


def test_engine_parsers_comet_init():
    input_file = Path(__file__).parent.parent / "data" / "BSA1_comet_2020_01_4.mzid"
    rt_lookup_path = (
        Path(__file__).parent.parent / "data" / "BSA1_ursgal_lookup.csv.bz2"
    )
    db_path = (
        Path(__file__).parent.parent / "data" / "test_Creinhardtii_target_decoy.fasta"
    )

    parser = Comet_2020_01_4_Parser(
        input_file,
        params={
            "rt_pickle_name": rt_lookup_path,
            "database": db_path,
            "modifications": [
                {
                    "aa": "M",
                    "type": "opt",
                    "position": "any",
                    "name": "Oxidation",
                },
                {
                    "aa": "C",
                    "type": "fix",
                    "position": "any",
                    "name": "Carbamidomethyl",
                },
                {
                    "aa": "*",
                    "type": "opt",
                    "position": "Prot-N-term",
                    "name": "Acetyl",
                },
            ],
            "omssa_mod_dir": Path(__file__).parent.parent / "data",
        },
    )


def test_engine_parsers_comet_file_matches_parser():
    msgf_parser_class = Comet_2020_01_4_Parser
    input_file = Path(__file__).parent.parent / "data" / "BSA1_comet_2020_01_4.mzid"
    assert msgf_parser_class.check_parser_compatibility(input_file) is True


def test_engine_parsers_comet_file_matches_parser_fail_with_omssa_file():
    msgf_parser_class = Comet_2020_01_4_Parser
    input_file = (
        Path(__file__).parent.parent
        / "data"
        / "test_Creinhardtii_QE_pH11_omssa_2_1_9.csv"
    )
    assert msgf_parser_class.check_parser_compatibility(input_file) is False


def test_engine_parsers_comet_unify():
    input_file = Path(__file__).parent.parent / "data" / "BSA1_comet_2020_01_4.mzid"
    rt_lookup_path = (
        Path(__file__).parent.parent / "data" / "BSA1_ursgal_lookup.csv.bz2"
    )
    db_path = (
        Path(__file__).parent.parent / "data" / "test_Creinhardtii_target_decoy.fasta"
    )

    parser = Comet_2020_01_4_Parser(
        input_file,
        params={
            "rt_pickle_name": rt_lookup_path,
            "database": db_path,
            "modifications": [
                {
                    "aa": "M",
                    "type": "opt",
                    "position": "any",
                    "name": "Oxidation",
                },
                {
                    "aa": "C",
                    "type": "fix",
                    "position": "any",
                    "name": "Carbamidomethyl",
                },
                {
                    "aa": "*",
                    "type": "opt",
                    "position": "Prot-N-term",
                    "name": "Acetyl",
                },
            ],
            "omssa_mod_dir": Path(__file__).parent.parent / "data",
        },
    )
    df = parser.unify()

    assert len(df) == 60


def test_engine_parsers_comet_get_peptide_lookup():
    input_file = Path(__file__).parent.parent / "data" / "BSA1_comet_2020_01_4.mzid"
    rt_lookup_path = (
        Path(__file__).parent.parent / "data" / "BSA1_ursgal_lookup.csv.bz2"
    )
    db_path = (
        Path(__file__).parent.parent / "data" / "test_Creinhardtii_target_decoy.fasta"
    )

    parser = Comet_2020_01_4_Parser(
        input_file,
        params={
            "rt_pickle_name": rt_lookup_path,
            "database": db_path,
            "modifications": [
                {
                    "aa": "M",
                    "type": "opt",
                    "position": "any",
                    "name": "Oxidation",
                },
                {
                    "aa": "C",
                    "type": "fix",
                    "position": "any",
                    "name": "Carbamidomethyl",
                },
                {
                    "aa": "*",
                    "type": "opt",
                    "position": "Prot-N-term",
                    "name": "Acetyl",
                },
            ],
            "omssa_mod_dir": Path(__file__).parent.parent / "data",
        },
    )
    lookup = parser._get_peptide_lookup()
    assert len(lookup) == 24
    assert "LVTDLTK;7:42.010565;" in lookup.keys()
    assert lookup["LVTDLTK;7:42.010565;"]["Sequence"] == "LVTDLTK"
    assert lookup["LVTDLTK;7:42.010565;"]["Modifications"] == "Acetyl:0"


def test_engine_parsers_comet_internal_next():
    input_file = Path(__file__).parent.parent / "data" / "BSA1_comet_2020_01_4.mzid"
    rt_lookup_path = (
        Path(__file__).parent.parent / "data" / "BSA1_ursgal_lookup.csv.bz2"
    )
    db_path = (
        Path(__file__).parent.parent / "data" / "test_Creinhardtii_target_decoy.fasta"
    )

    parser = Comet_2020_01_4_Parser(
        input_file,
        params={
            "rt_pickle_name": rt_lookup_path,
            "database": db_path,
            "modifications": [
                {
                    "aa": "M",
                    "type": "opt",
                    "position": "any",
                    "name": "Oxidation",
                },
                {
                    "aa": "C",
                    "type": "fix",
                    "position": "any",
                    "name": "Carbamidomethyl",
                },
                {
                    "aa": "*",
                    "type": "opt",
                    "position": "Prot-N-term",
                    "name": "Acetyl",
                },
            ],
            "omssa_mod_dir": Path(__file__).parent.parent / "data",
        },
    )
    df = parser.unify()
    row = df.iloc[1, :]
    assert (
        row["Raw data location"]
        == "/Users/cellzome/Dev/Gits/Ursgal/ursgal_master/example_data/BSA1.mzML"
    )
    assert row["Sequence"] == "LRCASIQK"
    # TODO: why arent fixed mods added
    assert row["Modifications"] == "Acetyl:0;Carbamidomethyl:3"
    assert row["Comet:spscore"] == "2.3000"
    assert row["Search Engine"] == "comet_2020_01_4"
