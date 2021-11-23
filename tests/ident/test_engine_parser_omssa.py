#!/usr/bin/env python
from pathlib import Path

import pytest

from unify_idents.engine_parsers.ident.omssa_2_1_9_parser import Omssa_Parser


def test_engine_parsers_omssa_init():
    input_file = pytest._test_path / "data" / "test_Creinhardtii_QE_pH11_omssa_2_1_9.csv"
    rt_lookup_path = pytest._test_path / "data" / "_ursgal_lookup.csv.bz2"
    db_path = pytest._test_path / "data" / "test_Creinhardtii_target_decoy.fasta"

    parser = Omssa_Parser(
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


def test_engine_parsers_omssa_check_parser_compatibility():
    input_file = pytest._test_path / "data" / "test_Creinhardtii_QE_pH11_omssa_2_1_9.csv"
    assert Omssa_Parser.check_parser_compatibility(input_file) is True


def test_engine_parsers_omssa_check_dataframe_integrity():
    input_file = pytest._test_path / "data" / "test_Creinhardtii_QE_pH11_omssa_2_1_9.csv"
    rt_lookup_path = pytest._test_path / "data" / "_ursgal_lookup.csv.bz2"
    db_path = pytest._test_path / "data" / "test_Creinhardtii_target_decoy.fasta"

    parser = Omssa_Parser(
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
            "Raw data location": "path/for/glory.mzML",
        },
    )
    df = parser.unify()
    assert pytest.approx(df["uCalc m/z"].mean()) == 826.67596
    assert (df["Raw data location"] == "path/for/glory.mzML").all()
    assert pytest.approx(df["Exp m/z"].mean()) == 826.7788

    assert df["Modifications"].str.contains("Acetyl:0").sum() == 0
    assert df["Modifications"].str.contains("Oxidation:").sum() == 93
    assert (
        df["Modifications"].str.count("Carbamidomethyl:")
        == df["Sequence"].str.count("C")
    ).all()
    assert df["Modifications"].str.count(":").sum() == 204


def test_replace_mod_strings():
    input_file = pytest._test_path / "data" / "test_Creinhardtii_QE_pH11_omssa_2_1_9.csv"
    rt_lookup_path = pytest._test_path / "data" / "_ursgal_lookup.csv.bz2"
    db_path = pytest._test_path / "data" / "test_Creinhardtii_target_decoy.fasta"

    parser = Omssa_Parser(
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
