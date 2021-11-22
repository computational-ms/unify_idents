#!/usr/bin/env python
from pathlib import Path
import pytest
from unify_idents.engine_parsers.ident.msgfplus_2021_03_22_parser import (
    MSGFPlus_2021_03_22_Parser,
)


def test_engine_parsers_msgfplus_init():
    input_file = Path(__file__).parent.parent / "data" / "BSA1_msgfplus_2021_03_22.mzid"
    rt_lookup_path = Path(__file__).parent.parent / "data" / "BSA1_ursgal_lookup.csv.bz2"
    db_path = (
        Path(__file__).parent.parent / "data" / "test_Creinhardtii_target_decoy.fasta"
    )

    parser = MSGFPlus_2021_03_22_Parser(
        input_file,
        params={
            "cpus": 2,
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


def test_engine_parsers_msgfplus_check_parser_compatibility():
    msgf_parser_class = MSGFPlus_2021_03_22_Parser
    input_file = Path(__file__).parent.parent / "data" / "BSA1_msgfplus_2021_03_22.mzid"
    assert msgf_parser_class.check_parser_compatibility(input_file) is True


def test_engine_parsers_msgfplus_check_parser_compatibility_fail_with_omsa_file():
    msgf_parser_class = MSGFPlus_2021_03_22_Parser
    input_file = (
        Path(__file__).parent.parent
        / "data"
        / "test_Creinhardtii_QE_pH11_omssa_2_1_9.csv"
    )
    assert msgf_parser_class.check_parser_compatibility(input_file) is False


def test_engine_parsers_msgfplus_check_dataframe_integrity():
    input_file = Path(__file__).parent.parent / "data" / "BSA1_msgfplus_2021_03_22.mzid"
    rt_lookup_path = Path(__file__).parent.parent / "data" / "BSA1_ursgal_lookup.csv.bz2"
    db_path = Path(__file__).parent.parent / "data" / "BSA.fasta"

    parser = MSGFPlus_2021_03_22_Parser(
        input_file,
        params={
            "cpus": 2,
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

    assert len(df) == 92
    assert pytest.approx(df["uCalc m/z"].mean()) == 488.03167


def test_engine_parsers_msgfplus_get_peptide_lookup():
    input_file = Path(__file__).parent.parent / "data" / "BSA1_msgfplus_2021_03_22.mzid"
    rt_lookup_path = Path(__file__).parent.parent / "data" / "BSA1_ursgal_lookup.csv.bz2"
    db_path = (
        Path(__file__).parent.parent / "data" / "test_Creinhardtii_target_decoy.fasta"
    )

    parser = MSGFPlus_2021_03_22_Parser(
        input_file,
        params={
            "cpus": 2,
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
    assert "Pep_YICDNQDTISSK" in lookup.keys()
    assert lookup["Pep_YICDNQDTISSK"]["Sequence"] == "YICDNQDTISSK"
    assert lookup["Pep_YICDNQDTISSK"]["Modifications"] == "Carbamidomethyl:3"


def test_get_single_spec_df():
    assert 1 == 2
