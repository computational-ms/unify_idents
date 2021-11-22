#!/usr/bin/env python
from pathlib import Path
import pytest
import pandas as pd

from unify_idents.engine_parsers.ident.xtandem_alanine import XTandemAlanine_Parser


def test_engine_parsers_xtandem_init():
    input_file = (
        Path(__file__).parent.parent
        / "data"
        / "test_Creinhardtii_QE_pH11_xtandem_alanine.xml"
    )
    rt_lookup_path = Path(__file__).parent.parent / "data" / "_ursgal_lookup.csv.bz2"
    db_path = (
        Path(__file__).parent.parent / "data" / "test_Creinhardtii_target_decoy.fasta"
    )

    parser = XTandemAlanine_Parser(
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
        },
    )


def test_engine_parsers_xtandem_file_matches_xtandem_parser():
    input_file = (
        Path(__file__).parent.parent
        / "data"
        / "test_Creinhardtii_QE_pH11_xtandem_alanine.xml"
    )

    assert XTandemAlanine_Parser.check_parser_compatibility(input_file) is True


def test_engine_parsers_msfragger_file_not_matches_xtandem_parser():
    input_file = (
        Path(__file__).parent.parent
        / "data"
        / "test_Creinhardtii_QE_pH11_mzml2mgf_0_0_1_msfragger_3.tsv"
    )

    assert XTandemAlanine_Parser.check_parser_compatibility(input_file) is False


def test_engine_parsers_xtandem_check_dataframe_integrity():
    input_file = (
        Path(__file__).parent.parent
        / "data"
        / "test_Creinhardtii_QE_pH11_xtandem_alanine.xml"
    )
    rt_lookup_path = Path(__file__).parent.parent / "data" / "_ursgal_lookup.csv.bz2"
    db_path = (
        Path(__file__).parent.parent / "data" / "test_Creinhardtii_target_decoy.fasta"
    )

    parser = XTandemAlanine_Parser(
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
        },
    )
    df = parser.unify()
    assert len(parser.root) == 79
    assert pytest.approx(df["uCalc m/z"].mean(), 457.85944)


def test_get_single_spec_df():
    assert 1 == 2


def test_map_mod_names():
    assert 1 == 2


def test_engine_parsers_xtandem_nterminal_mod():
    input_file = (
        Path(__file__).parent.parent
        / "data"
        / "test_Creinhardtii_QE_pH11_xtandem_alanine.xml"
    )
    rt_lookup_path = Path(__file__).parent.parent / "data" / "_ursgal_lookup.csv.bz2"
    db_path = (
        Path(__file__).parent.parent / "data" / "test_Creinhardtii_target_decoy.fasta"
    )

    parser = XTandemAlanine_Parser(
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
            "Raw file location": "test_Creinhardtii_QE_pH11.mzML",
            "15N": False,
        },
    )
    df = parser.unify()
    relevant_row = df[df["Sequence"] == "WGLVSSELQTSEAETPGLK"]
    assert relevant_row["Modifications"].tolist() == ["Acetyl:0"]


def test_engine_parsers_xtandem_multiple_psms():
    input_file = Path(__file__).parent.parent / "data" / "multiple_psms_xtandem.xml"
    rt_lookup_path = Path(__file__).parent.parent / "data" / "_ursgal_lookup.csv.bz2"
    db_path = (
        Path(__file__).parent.parent / "data" / "human_ecoli_test_target_decoy.fasta"
    )

    parser = XTandemAlanine_Parser(
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
            "Raw file location": "test_Creinhardtii_QE_pH11.mzML",
            "15N": False,
        },
    )

    # Test file with:
    #    - one sequence in first group
    #    - 3 sequences in second group

    df = parser.unify()
    assert len(df) == 4

    assert set(df["Sequence"]) == {
        "ITIPITLRMLIAK",
        "SMMNGGSSPESDVGTDNK",
        "SMMNGGSSPESDVGTDNK",
        "SMMNGGSSPESDVGTDNK",
    }
    assert set(df["Spectrum ID"]) == {12833, 14525}
    assert set(df["Modifications"]) == {
        "Acetyl:0",
        "",
        "Oxidation:2",
        "Oxidation:3",
    }


def test_engine_parsers_xtandem_map_mod_names():
    input_file = (
        Path(__file__).parent.parent
        / "data"
        / "test_Creinhardtii_QE_pH11_xtandem_alanine.xml"
    )
    rt_lookup_path = Path(__file__).parent.parent / "data" / "_ursgal_lookup.csv.bz2"
    db_path = (
        Path(__file__).parent.parent / "data" / "test_Creinhardtii_target_decoy.fasta"
    )

    parser = XTandemAlanine_Parser(
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
        },
    )
    test_df = pd.DataFrame({"Modifications": [["57.021464:0"]], "Sequence": ["CERK"]})
    assert parser.map_mod_names(test_df)["Modifications"][0] == "Carbamidomethyl:1"


def test_engine_parsers_xtandem_map_mod_names_nterm():
    input_file = (
        Path(__file__).parent.parent
        / "data"
        / "test_Creinhardtii_QE_pH11_xtandem_alanine.xml"
    )
    rt_lookup_path = Path(__file__).parent.parent / "data" / "_ursgal_lookup.csv.bz2"
    db_path = (
        Path(__file__).parent.parent / "data" / "test_Creinhardtii_target_decoy.fasta"
    )

    parser = XTandemAlanine_Parser(
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
        },
    )

    row = pd.DataFrame(
        {"Modifications": [["57.021464:0", "42.010565:0"]], "Sequence": ["CERK"]}
    )
    assert set(parser.map_mod_names(row)["Modifications"][0].split(";")) == {
        "Carbamidomethyl:1",
        "Acetyl:0",
    }
