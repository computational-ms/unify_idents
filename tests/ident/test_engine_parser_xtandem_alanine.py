#!/usr/bin/env python
from pathlib import Path

import pandas as pd
import pytest

import xml.etree.ElementTree as ETree
from unify_idents.engine_parsers.ident.xtandem_alanine import (
    XTandemAlanine_Parser,
    _get_single_spec_df,
)


def test_engine_parsers_xtandem_init():
    input_file = (
        pytest._test_path / "data" / "test_Creinhardtii_QE_pH11_xtandem_alanine.xml"
    )
    rt_lookup_path = pytest._test_path / "data" / "_ursgal_lookup.csv.bz2"
    db_path = pytest._test_path / "data" / "test_Creinhardtii_target_decoy.fasta"

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
        pytest._test_path / "data" / "test_Creinhardtii_QE_pH11_xtandem_alanine.xml"
    )

    assert XTandemAlanine_Parser.check_parser_compatibility(input_file) is True


def test_engine_parsers_msfragger_file_not_matches_xtandem_parser():
    input_file = pytest._test_path / "data" / "test_Creinhardtii_QE_pH11_msfragger_3.tsv"

    assert XTandemAlanine_Parser.check_parser_compatibility(input_file) is False


def test_engine_parsers_xtandem_check_dataframe_integrity():
    input_file = (
        pytest._test_path / "data" / "test_Creinhardtii_QE_pH11_xtandem_alanine.xml"
    )
    rt_lookup_path = pytest._test_path / "data" / "_ursgal_lookup.csv.bz2"
    db_path = pytest._test_path / "data" / "test_Creinhardtii_target_decoy.fasta"

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
    assert (df["Raw data location"] == "path/for/glory.mzML").all()
    assert pytest.approx(df["uCalc m/z"].mean()) == 796.4324
    assert pytest.approx(df["Exp m/z"].mean()) == 796.71967

    assert df["Modifications"].str.contains("Acetyl:0").sum() == 1
    assert df["Modifications"].str.contains("Oxidation:").sum() == 23
    assert (
        df["Modifications"].str.count("Carbamidomethyl:")
        == df["Sequence"].str.count("C")
    ).all()
    assert df["Modifications"].str.count(":").sum() == 50


def test_get_single_spec_df():
    input_file = (
        pytest._test_path / "data" / "test_Creinhardtii_QE_pH11_xtandem_alanine.xml"
    )
    element = ETree.parse(input_file).getroot()[0]
    ref_dict = {
        "Exp m/z": None,
        "Calc m/z": None,
        "Spectrum Title": None,
        "Raw data location": "path/for/glory.mgf",
        "Search Engine": "xtandem_alanine",
        "Spectrum ID": None,
        "Modifications": None,
        "Retention Time (s)": None,
        "X!Tandem:delta": None,
        "X!Tandem:nextscore": None,
        "X!Tandem:y_score": None,
        "X!Tandem:y_ions": None,
        "X!Tandem:b_score": None,
        "X!Tandem:b_ions": None,
        "Sequence": None,
        "Charge": None,
        "X!Tandem:Hyperscore": None,
    }
    mapping_dict = {
        "delta": "X!Tandem:delta",
        "nextscore": "X!Tandem:nextscore",
        "y_score": "X!Tandem:y_score",
        "y_ions": "X!Tandem:y_ions",
        "b_score": "X!Tandem:b_score",
        "b_ions": "X!Tandem:b_ions",
        "seq": "Sequence",
        "z": "Charge",
        "hyperscore": "X!Tandem:Hyperscore",
    }

    result = _get_single_spec_df(ref_dict, mapping_dict, element)

    assert isinstance(result, pd.DataFrame)
    assert (
        result.values
        == [
            None,
            "1315.5700",
            "test_Creinhardtii_QE_pH11.10381.10381.3",
            "path/for/glory.mgf",
            "xtandem_alanine",
            "10381",
            list(["15.99492:5"]),
            None,
            "0.0057",
            "8.0",
            "9.9",
            "5",
            "0.0",
            "0",
            "DDVHNMGADGIR",
            "3",
            "14.2",
        ]
    ).all()


def test_engine_parsers_xtandem_nterminal_mod():
    input_file = (
        pytest._test_path / "data" / "test_Creinhardtii_QE_pH11_xtandem_alanine.xml"
    )
    rt_lookup_path = pytest._test_path / "data" / "_ursgal_lookup.csv.bz2"
    db_path = pytest._test_path / "data" / "test_Creinhardtii_target_decoy.fasta"

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
    input_file = pytest._test_path / "data" / "multiple_psms_xtandem.xml"
    rt_lookup_path = pytest._test_path / "data" / "_ursgal_lookup.csv.bz2"
    db_path = pytest._test_path / "data" / "human_ecoli_test_target_decoy.fasta"

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
        pytest._test_path / "data" / "test_Creinhardtii_QE_pH11_xtandem_alanine.xml"
    )
    rt_lookup_path = pytest._test_path / "data" / "_ursgal_lookup.csv.bz2"
    db_path = pytest._test_path / "data" / "test_Creinhardtii_target_decoy.fasta"

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
        pytest._test_path / "data" / "test_Creinhardtii_QE_pH11_xtandem_alanine.xml"
    )
    rt_lookup_path = pytest._test_path / "data" / "_ursgal_lookup.csv.bz2"
    db_path = pytest._test_path / "data" / "test_Creinhardtii_target_decoy.fasta"

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
