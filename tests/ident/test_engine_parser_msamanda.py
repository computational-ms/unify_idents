#!/usr/bin/env python

import pytest

from unify_idents.engine_parsers.ident.msamanda_2_parser import MSAmanda_2_Parser


def test_engine_parsers_msamanda_init():
    input_file = pytest._test_path / "data" / "BSA1_msamanda_2_0_0_17442.csv"

    parser = MSAmanda_2_Parser(
        input_file,
        params={
            "cpus": 2,
            "enzyme": "(?<=[KR])(?![P])",
            "terminal_cleavage_site_integrity": "any",
            "validation_score_field": {"msamanda_2_0_0_17442": "amanda:score"},
            "bigger_scores_better": {"msamanda_2_0_0_17442": True},
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
        },
    )


def test_engine_parsers_msamanda_check_parser_compatibility():
    input_file = pytest._test_path / "data" / "BSA1_msamanda_2_0_0_17442.csv"
    assert MSAmanda_2_Parser.check_parser_compatibility(input_file) is True


def test_engine_parsers_msamanda_check_dataframe_integrity():
    input_file = pytest._test_path / "data" / "BSA1_msamanda_2_0_0_17442.csv"
    rt_lookup_path = pytest._test_path / "data" / "BSA1_ursgal_lookup.csv"
    db_path = pytest._test_path / "data" / "BSA.fasta"

    parser = MSAmanda_2_Parser(
        input_file,
        params={
            "cpus": 2,
            "rt_pickle_name": rt_lookup_path,
            "database": db_path,
            "enzyme": "(?<=[KR])(?![P])",
            "terminal_cleavage_site_integrity": "any",
            "validation_score_field": {"msamanda_2_0_0_17442": "amanda:score"},
            "bigger_scores_better": {"msamanda_2_0_0_17442": True},
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
        },
    )
    df = parser.unify()
    assert len(df) == 87
    assert pytest.approx(df["ucalc_mz"].mean()) == 485.2679
    assert pytest.approx(df["exp_mz"].mean()) == 485.26797

    assert df["modifications"].str.contains("Acetyl:0").sum() == 0
    assert (
        df["modifications"].str.count("Carbamidomethyl:")
        == df["sequence"].str.count("C")
    ).all()
    assert df["modifications"].str.count(":").sum() == 66
    assert (df["raw_data_location"] == "path/for/glory.mzML").all()


def test_map_mod_translation():
    input_file = pytest._test_path / "data" / "BSA1_msamanda_2_0_0_17442.csv"

    parser = MSAmanda_2_Parser(
        input_file,
        params={
            "cpus": 2,
            "enzyme": "(?<=[KR])(?![P])",
            "terminal_cleavage_site_integrity": "any",
            "validation_score_field": {"msamanda_2_0_0_17442": "amanda:score"},
            "bigger_scores_better": {"msamanda_2_0_0_17442": True},
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
        },
    )
    converted = parser._map_mod_translation(row=["C3(Carbamidomethyl|57.021464|fixed)"])
    assert converted == "Carbamidomethyl:3;"
