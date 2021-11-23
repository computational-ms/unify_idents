#!/usr/bin/env python
from pathlib import Path

import pandas as pd
import pytest
import xml.etree.ElementTree as ETree

from unify_idents.engine_parsers.ident.comet_2020_01_4_parser import (
    Comet_2020_01_4_Parser,
    _get_single_spec_df,
)


def test_engine_parsers_comet_init():
    input_file = pytest._test_path / "data" / "BSA1_comet_2020_01_4.mzid"
    rt_lookup_path = pytest._test_path / "data" / "BSA1_ursgal_lookup.csv.bz2"
    db_path = pytest._test_path / "data" / "test_Creinhardtii_target_decoy.fasta"

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
        },
    )


def test_engine_parsers_comet_check_parser_compatibility():
    msgf_parser_class = Comet_2020_01_4_Parser
    input_file = pytest._test_path / "data" / "BSA1_comet_2020_01_4.mzid"
    assert msgf_parser_class.check_parser_compatibility(input_file) is True


def test_engine_parsers_comet_check_parser_compatibility_fail_with_omssa_file():
    msgf_parser_class = Comet_2020_01_4_Parser
    input_file = pytest._test_path / "data" / "test_Creinhardtii_QE_pH11_omssa_2_1_9.csv"
    assert msgf_parser_class.check_parser_compatibility(input_file) is False


def test_engine_parsers_comet_check_dataframe_integrity():
    input_file = pytest._test_path / "data" / "BSA1_comet_2020_01_4.mzid"
    rt_lookup_path = pytest._test_path / "data" / "BSA1_ursgal_lookup.csv.bz2"
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
        },
    )
    df = parser.unify()
    assert pytest.approx(df["uCalc m/z"].mean()) == 457.85944
    assert pytest.approx(df["Exp m/z"].mean()) == 457.87625

    assert df["Modifications"].str.contains("Acetyl:0").sum() == 5
    assert df["Modifications"].str.contains("Oxidation:").sum() == 0
    assert (
        df["Modifications"].str.count("Carbamidomethyl:")
        == df["Sequence"].str.count("C")
    ).all()
    assert df["Modifications"].str.count(":").sum() == 38
    assert (df["Raw data location"] == "path/for/glory.mzML").all()


def test_get_single_spec_df():
    input_file = pytest._test_path / "data" / "BSA1_comet_2020_01_4.mzid"
    element = (
        ETree.parse(input_file)
        .getroot()
        .find(".//{*}SpectrumIdentificationList/{*}SpectrumIdentificationResult")
    )
    ref_dict = {
        "Exp m/z": None,
        "Calc m/z": None,
        "Spectrum Title": None,
        "Raw data location": "path/for/glory.mgf",
        "Search Engine": "comet_2020_01_4",
        "Spectrum ID": None,
        "Modifications": None,
        "Retention Time (s)": None,
        "Charge": None,
        "Comet:Score": None,
        "Comet:DeltaCn": None,
        "Comet:XCorr": None,
        "Comet:EValue": None,
        "Sequence": None,
        "Comet:SpecEValue": None,
        "Comet:Num Matched Ions": None,
        "Comet:Num Unmatched Ions": None,
    }
    mapping_dict = {
        "chargeState": "Charge",
        "Comet:spscore": "Comet:Score",
        "Comet:deltacn": "Comet:DeltaCn",
        "Comet:xcorr": "Comet:XCorr",
        "Comet:expectation value": "Comet:EValue",
        "peptide_ref": "Sequence",
        "experimentalMassToCharge": "Exp m/z",
        "calculatedMassToCharge": "Calc m/z",
        "SpecEValue": "Comet:SpecEValue",
        "number of matched peaks": "Comet:Num Matched Ions",
        "number of unmatched peaks": "Comet:Num Unmatched Ions",
    }

    result = _get_single_spec_df(ref_dict, mapping_dict, element)

    assert isinstance(result, pd.DataFrame)
    assert (
        result.values
        == [
            [
                "358.174682",
                "358.174575",
                None,
                "path/for/glory.mgf",
                "comet_2020_01_4",
                "2458",
                None,
                None,
                "3",
                "5.9000",
                "1.0000",
                "0.2825",
                "3.76E+01",
                "SHCIAEVEK;",
                None,
                "3",
                "29",
            ]
        ]
    ).all()
