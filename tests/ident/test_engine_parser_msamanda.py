#!/usr/bin/env python
from pathlib import Path

import pytest

from unify_idents import Unify
from unify_idents.engine_parsers.ident.msamanda_2_parser import MSAmanda_2_Parser


def test_engine_parsers_msamanda_init():
    input_file = Path(__file__).parent.parent / "data" / "BSA1_msamanda_2_0_0_17442.csv"
    rt_lookup_path = (
        Path(__file__).parent.parent / "data" / "BSA1_ursgal_lookup.csv.bz2"
    )
    db_path = Path(__file__).parent / "data" / "BSA.fasta"

    parser = MSAmanda_2_Parser(
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
            # "msamanda_mod_dir": Path(__file__).parent / "data",
        },
    )


def test_engine_parsers_msamanda_check_parser_compatibility():
    input_file = Path(__file__).parent.parent / "data" / "BSA1_msamanda_2_0_0_17442.csv"
    assert MSAmanda_2_Parser.check_parser_compatibility(input_file) is True


def test_engine_parsers_msamanda_iterable():
    input_file = Path(__file__).parent.parent / "data" / "BSA1_msamanda_2_0_0_17442.csv"
    rt_lookup_path = (
        Path(__file__).parent.parent / "data" / "BSA1_ursgal_lookup.csv.bz2"
    )
    db_path = Path(__file__).parent.parent / "data" / "BSA.fasta"

    parser = MSAmanda_2_Parser(
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
            "Raw data location": "/Users/cellzome/Dev/Gits/Ursgal/ursgal_master/example_data/test_Creinhardtii_QE_pH11.mzML",
        },
    )
    df = parser.unify()
    assert len(df) == 87


def test_engine_parsers_msamanda_next():
    input_file = Path(__file__).parent.parent / "data" / "BSA1_msamanda_2_0_0_17442.csv"
    rt_lookup_path = (
        Path(__file__).parent.parent / "data" / "BSA1_ursgal_lookup.csv.bz2"
    )
    db_path = Path(__file__).parent.parent / "data" / "BSA.fasta"

    parser = Unify(
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
            "omssa_mod_dir": Path(__file__).parent / "data",
            "Raw data location": "/Users/cellzome/Dev/Gits/Ursgal/ursgal_master/example_data/test_Creinhardtii_QE_pH11.mzML",
        },
    )
    df = parser.get_dataframe()
    first_row = df.iloc[0, :]
    assert first_row["Sequence"] == "SHCIAEVEK"
    assert first_row["Modifications"] == "Carbamidomethyl:3"
    assert (
        first_row["Raw data location"]
        == "/Users/cellzome/Dev/Gits/Ursgal/ursgal_master/example_data/test_Creinhardtii_QE_pH11.mzML"
    )
    assert first_row["Charge"] == 3
    assert float(first_row["Amanda:Score"]) == pytest.approx(83.16219696000763)
    assert float(first_row["Exp m/z"]) == pytest.approx(358.174682617188)
