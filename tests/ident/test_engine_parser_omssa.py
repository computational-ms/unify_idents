#!/usr/bin/env python

import pytest

from unify_idents.engine_parsers.ident.omssa_2_1_9_parser import Omssa_Parser


def test_engine_parsers_omssa_init():
    input_file = (
        pytest._test_path / "data" / "test_Creinhardtii_QE_pH11_omssa_2_1_9.csv"
    )

    parser = Omssa_Parser(
        input_file,
        params={
            "cpus": 2,
            "enzyme": "(?<=[KR])(?![P])",
            "terminal_cleavage_site_integrity": "any",
            "validation_score_field": {"omssa_2_1_9": "omssa:pvalue"},
            "bigger_scores_better": {"omssa_2_1_9": False},
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
            "xml_file_list": [
                pytest._test_path / "data" / "mods.xml",
                pytest._test_path / "data" / "usermods.xml",
            ],
        },
    )


def test_engine_parsers_omssa_check_parser_compatibility():
    input_file = (
        pytest._test_path / "data" / "test_Creinhardtii_QE_pH11_omssa_2_1_9.csv"
    )
    assert Omssa_Parser.check_parser_compatibility(input_file) is True


def test_engine_parsers_omssa_check_dataframe_integrity():
    input_file = (
        pytest._test_path / "data" / "test_Creinhardtii_QE_pH11_omssa_2_1_9.csv"
    )
    rt_lookup_path = pytest._test_path / "data" / "_ursgal_lookup.csv"
    db_path = pytest._test_path / "data" / "test_Creinhardtii_target_decoy.fasta"

    parser = Omssa_Parser(
        input_file,
        params={
            "cpus": 2,
            "rt_pickle_name": rt_lookup_path,
            "database": db_path,
            "enzyme": "(?<=[KR])(?![P])",
            "terminal_cleavage_site_integrity": "any",
            "validation_score_field": {"omssa_2_1_9": "omssa:pvalue"},
            "bigger_scores_better": {"omssa_2_1_9": False},
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
            "xml_file_list": [
                pytest._test_path / "data" / "mods.xml",
                pytest._test_path / "data" / "usermods.xml",
            ],
        },
    )
    df = parser.unify()
    assert pytest.approx(df["ucalc_mz"].mean()) == 826.7793
    assert (df["raw_data_location"] == "path/for/glory.mzML").all()
    assert pytest.approx(df["exp_mz"].mean()) == 826.7788

    assert df["modifications"].str.contains("Acetyl:0").sum() == 5
    assert df["modifications"].str.contains("Oxidation:").sum() == 93
    assert (
        df["modifications"].str.count("Carbamidomethyl:")
        == df["sequence"].str.count("C")
    ).all()
    assert df["modifications"].str.count(":").sum() == 204


def test_replace_mod_strings():
    input_file = (
        pytest._test_path / "data" / "test_Creinhardtii_QE_pH11_omssa_2_1_9.csv"
    )

    parser = Omssa_Parser(
        input_file,
        params={
            "cpus": 2,
            "enzyme": "(?<=[KR])(?![P])",
            "terminal_cleavage_site_integrity": "any",
            "validation_score_field": {"omssa_2_1_9": "omssa:pvalue"},
            "bigger_scores_better": {"omssa_2_1_9": False},
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
            "xml_file_list": [
                pytest._test_path / "data" / "mods.xml",
                pytest._test_path / "data" / "usermods.xml",
            ],
        },
    )
    raw_mod = "oxidation of M:4"
    converted = parser._replace_mod_strings(
        raw_mod,
        {
            "acetylation of protein n-term": {
                "unimod_id": None,
                "omssa_unimod_id": "1",
                "unimod_name": "Acetyl",
                "omssa_name": "acetylation of protein n-term",
                "aa_targets": [],
            },
            "oxidation of M": {
                "unimod_id": None,
                "omssa_unimod_id": "35",
                "unimod_name": "Oxidation",
                "omssa_name": "oxidation of M",
                "aa_targets": ["M"],
            },
        },
    )
    assert converted == "Oxidation:4"
