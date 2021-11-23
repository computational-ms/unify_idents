#!/usr/bin/env python
from pathlib import Path

import xml.etree.ElementTree as ETree
import pandas as pd
import pytest

from unify_idents.engine_parsers.ident.msgfplus_2021_03_22_parser import (
    MSGFPlus_2021_03_22_Parser,
    _get_single_spec_df,
)


def test_engine_parsers_msgfplus_init():
    input_file = pytest._test_path / "data" / "BSA1_msgfplus_2021_03_22.mzid"
    rt_lookup_path = pytest._test_path / "data" / "BSA1_ursgal_lookup.csv.bz2"
    db_path = pytest._test_path / "data" / "test_Creinhardtii_target_decoy.fasta"

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
            "omssa_mod_dir": pytest._test_path / "data",
        },
    )


def test_engine_parsers_msgfplus_check_parser_compatibility():
    msgf_parser_class = MSGFPlus_2021_03_22_Parser
    input_file = pytest._test_path / "data" / "BSA1_msgfplus_2021_03_22.mzid"
    assert msgf_parser_class.check_parser_compatibility(input_file) is True


def test_engine_parsers_msgfplus_check_parser_compatibility_fail_with_omssa_file():
    msgf_parser_class = MSGFPlus_2021_03_22_Parser
    input_file = pytest._test_path / "data" / "test_Creinhardtii_QE_pH11_omssa_2_1_9.csv"
    assert msgf_parser_class.check_parser_compatibility(input_file) is False


def test_engine_parsers_msgfplus_check_dataframe_integrity():
    input_file = pytest._test_path / "data" / "BSA1_msgfplus_2021_03_22.mzid"
    rt_lookup_path = pytest._test_path / "data" / "BSA1_ursgal_lookup.csv.bz2"
    db_path = pytest._test_path / "data" / "BSA.fasta"

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
            "omssa_mod_dir": pytest._test_path / "data",
        },
    )
    df = parser.unify()
    assert pytest.approx(df["Exp m/z"].mean()) == 488.0319
    assert len(df) == 92
    assert pytest.approx(df["uCalc m/z"].mean()) == 488.03167
    assert (df["Raw data location"] == "path/for/glory.mzML").all()
    assert df["Modifications"].str.contains("Acetyl:0").sum() == 0
    assert df["Modifications"].str.contains("Oxidation:").sum() == 0
    assert (
        df["Modifications"].str.count("Carbamidomethyl:")
        == df["Sequence"].str.count("C")
    ).all()
    assert df["Modifications"].str.count(":").sum() == 71


def test_engine_parsers_msgfplus_get_peptide_lookup():
    input_file = pytest._test_path / "data" / "BSA1_msgfplus_2021_03_22.mzid"
    rt_lookup_path = pytest._test_path / "data" / "BSA1_ursgal_lookup.csv.bz2"
    db_path = pytest._test_path / "data" / "test_Creinhardtii_target_decoy.fasta"

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
            "omssa_mod_dir": pytest._test_path / "data",
        },
    )
    lookup = parser._get_peptide_lookup()
    assert len(lookup) == 24
    assert "Pep_YICDNQDTISSK" in lookup.keys()
    assert lookup["Pep_YICDNQDTISSK"]["Sequence"] == "YICDNQDTISSK"
    assert lookup["Pep_YICDNQDTISSK"]["Modifications"] == "Carbamidomethyl:3"


def test_get_single_spec_df():
    input_file = pytest._test_path / "data" / "BSA1_msgfplus_2021_03_22.mzid"
    element = (
        ETree.parse(input_file).getroot().find(".//{*}SpectrumIdentificationResult")
    )
    ref_dict = {
        "Exp m/z": None,
        "Calc m/z": None,
        "Spectrum Title": None,
        "Raw data location": "path/for/glory.mzML",
        "Search Engine": "msgfplus_2021_03_22",
        "Spectrum ID": None,
        "Modifications": None,
        "Retention Time (s)": None,
        "Charge": None,
        "MS-GF:DeNovoScore": None,
        "MS-GF:EValue": None,
        "MS-GF:RawScore": None,
        "Sequence": None,
        "MS-GF:SpecEValue": None,
        "MS-GF:Num Matched Ions": None,
    }
    mapping_dict = {
        "chargeState": "Charge",
        "MS-GF:DeNovoScore": "MS-GF:DeNovoScore",
        "MS-GF:EValue": "MS-GF:EValue",
        "MS-GF:RawScore": "MS-GF:RawScore",
        "peptide_ref": "Sequence",
        "experimentalMassToCharge": "Exp m/z",
        "calculatedMassToCharge": "Calc m/z",
        "scan number(s)": "Spectrum ID",
        "MS-GF:SpecEValue": "MS-GF:SpecEValue",
        "spectrum title": "Spectrum Title",
        "NumMatchedMainIons": "MS-GF:Num Matched Ions",
    }

    result = _get_single_spec_df(ref_dict, mapping_dict, element)

    assert isinstance(result, pd.DataFrame)
    assert (
        result.values
        == [
            [
                "722.3272094726562",
                "722.3246459960938",
                "glory.2791.2791.2",
                "path/for/glory.mzML",
                "msgfplus_2021_03_22",
                "2791",
                None,
                None,
                "2",
                "40",
                "2.6986221E-12",
                "40",
                "Pep_YICDNQDTISSK",
                "4.4458354E-15",
                "3",
            ]
        ]
    ).all()
