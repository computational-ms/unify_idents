#!/usr/bin/env python
from pathlib import Path

import pytest

from unify_idents.engine_parsers.ident.comet_2020_01_4_parser import (
    Comet_2020_01_4_Parser,
)


def test_engine_parsers_comet_init():
    input_file = pytest._test_path / "data" / "BSA1_comet_2020_01_4.mzid"
    rt_lookup_path = (
        pytest._test_path / "data" / "BSA1_ursgal_lookup.csv.bz2"
    )
    db_path = (
        pytest._test_path / "data" / "test_Creinhardtii_target_decoy.fasta"
    )

    parser = Comet_2020_01_4_Parser(
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
            "omssa_mod_dir": pytest._test_path / "data",
        },
    )


def test_engine_parsers_comet_check_parser_compatibility():
    msgf_parser_class = Comet_2020_01_4_Parser
    input_file = pytest._test_path / "data" / "BSA1_comet_2020_01_4.mzid"
    assert msgf_parser_class.check_parser_compatibility(input_file) is True


def test_engine_parsers_comet_check_parser_compatibility_fail_with_omssa_file():
    msgf_parser_class = Comet_2020_01_4_Parser
    input_file = (
        pytest._test_path
        / "data"
        / "test_Creinhardtii_QE_pH11_omssa_2_1_9.csv"
    )
    assert msgf_parser_class.check_parser_compatibility(input_file) is False


def test_engine_parsers_comet_check_dataframe_integrity():
    input_file = pytest._test_path / "data" / "BSA1_comet_2020_01_4.mzid"
    rt_lookup_path = (
        pytest._test_path / "data" / "BSA1_ursgal_lookup.csv.bz2"
    )
    db_path = pytest._test_path / "data" / "BSA.fasta"

    parser = Comet_2020_01_4_Parser(
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
            "omssa_mod_dir": pytest._test_path / "data",
        },
    )
    df = parser.unify()

    assert len(df) == 60
    assert pytest.approx(df["uCalc m/z"].mean()) == 485.26791

    assert df["Modifications"].str.contains("Acetyl:0").sum() == 5
    assert (df["Modifications"].str.count("Carbamidomethyl:") == df["Sequence"].str.count("C")).all()
    assert df["Modifications"].str.count(":").sum() == 38

    # assert mean uCalc mz
    # assert mean Exp mz


def test_get_single_spec_df():
    assert 1 == 2


def test_map_mods_and_sequences():
    assert 1 == 2
