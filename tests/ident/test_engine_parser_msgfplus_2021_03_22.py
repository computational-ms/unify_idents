#!/usr/bin/env python
import pandas as pd
import pytest
from lxml import etree

from unify_idents.engine_parsers.ident.msgfplus_2021_03_22_parser import (
    MSGFPlus_2021_03_22_Parser,
    _get_single_spec_df,
)


def test_engine_parsers_msgfplus_init():
    input_file = pytest._test_path / "data" / "BSA1_msgfplus_2021_03_22.mzid"
    parser = MSGFPlus_2021_03_22_Parser(
        input_file,
        params={
            "cpus": 2,
            "enzyme": "(?<=[KR])(?![P])",
            "terminal_cleavage_site_integrity": "any",
            "validation_score_field": {"msgfplus_2021_03_22": "ms-gf:spec_evalue"},
            "bigger_scores_better": {"msgfplus_2021_03_22": False},
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


def test_engine_parsers_msgfplus_check_parser_compatibility():
    msgf_parser_class = MSGFPlus_2021_03_22_Parser
    input_file = pytest._test_path / "data" / "BSA1_msgfplus_2021_03_22.mzid"
    assert msgf_parser_class.check_parser_compatibility(input_file) is True


def test_engine_parsers_msgfplus_check_parser_compatibility_fail_with_omssa_file():
    msgf_parser_class = MSGFPlus_2021_03_22_Parser
    input_file = (
        pytest._test_path / "data" / "test_Creinhardtii_QE_pH11_omssa_2_1_9.csv"
    )
    assert msgf_parser_class.check_parser_compatibility(input_file) is False


def test_engine_parsers_msgfplus_check_dataframe_integrity():
    input_file = pytest._test_path / "data" / "BSA1_msgfplus_2021_03_22.mzid"
    rt_lookup_path = pytest._test_path / "data" / "BSA1_ursgal_lookup.csv"
    db_path = pytest._test_path / "data" / "BSA.fasta"

    parser = MSGFPlus_2021_03_22_Parser(
        input_file,
        params={
            "cpus": 2,
            "rt_pickle_name": rt_lookup_path,
            "database": db_path,
            "enzyme": "(?<=[KR])(?![P])",
            "terminal_cleavage_site_integrity": "any",
            "validation_score_field": {"msgfplus_2021_03_22": "ms-gf:spec_evalue"},
            "bigger_scores_better": {"msgfplus_2021_03_22": False},
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
    print(df)
    assert pytest.approx(df["exp_mz"].mean()) == 488.0319
    assert len(df) == 92
    assert pytest.approx(df["ucalc_mz"].mean()) == 488.0288
    assert (df["raw_data_location"] == "path/for/glory.mzML").all()
    assert df["modifications"].str.contains("Acetyl:0").sum() == 0
    assert df["modifications"].str.contains("Oxidation:").sum() == 0
    assert (
        df["modifications"].str.count("Carbamidomethyl:")
        == df["sequence"].str.count("C")
    ).all()
    assert df["modifications"].str.count(":").sum() == 71


def test_engine_parsers_msgfplus_get_peptide_lookup():
    input_file = pytest._test_path / "data" / "BSA1_msgfplus_2021_03_22.mzid"

    parser = MSGFPlus_2021_03_22_Parser(
        input_file,
        params={
            "cpus": 2,
            "enzyme": "(?<=[KR])(?![P])",
            "terminal_cleavage_site_integrity": "any",
            "validation_score_field": {"msgfplus_2021_03_22": "ms-gf:spec_evalue"},
            "bigger_scores_better": {"msgfplus_2021_03_22": False},
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
    lookup = parser._get_peptide_lookup()
    assert len(lookup) == 24
    assert "Pep_YICDNQDTISSK" in lookup.keys()
    assert lookup["Pep_YICDNQDTISSK"]["sequence"] == "YICDNQDTISSK"
    assert lookup["Pep_YICDNQDTISSK"]["modifications"] == "Carbamidomethyl:3"


def test_get_single_spec_df():
    input_file = pytest._test_path / "data" / "BSA1_msgfplus_2021_03_22.mzid"
    element = (
        etree.parse(input_file).getroot().find(".//{*}SpectrumIdentificationResult")
    )
    ref_dict = {
        "exp_mz": None,
        "calc_mz": None,
        "spectrum_title": None,
        "search_engine": "msgfplus_2021_03_22",
        "spectrum_id": None,
        "modifications": None,
        "retention time_seconds": None,
        "charge": None,
        "ms-gf:denovoscore": None,
        "ms-gf:evalue": None,
        "ms-gf:raw_score": None,
        "sequence": None,
        "ms-gf:spec_evalue": None,
        "ms-gf:num_matched_ions": None,
    }
    mapping_dict = {
        "chargeState": "charge",
        "MS-GF:DeNovoScore": "ms-gf:denovoscore",
        "MS-GF:EValue": "ms-gf:evalue",
        "MS-GF:RawScore": "ms-gf:rawscore",
        "peptide_ref": "sequence",
        "experimentalMassToCharge": "exp_mz",
        "calculatedMassToCharge": "calc_mz",
        "scan number(s)": "spectrum_id",
        "MS-GF:SpecEValue": "ms-gf:spec_evalue",
        "spectrum title": "spectrum_title",
        "NumMatchedMainIons": "ms-gf:num_matched_ions",
    }
    _get_single_spec_df.reference_dict = ref_dict
    _get_single_spec_df.mapping_dict = mapping_dict
    result = _get_single_spec_df(etree.tostring(element))

    assert isinstance(result, pd.DataFrame)
    assert (
        result.values
        == [
            [
                "722.3272094726562",
                "722.3246459960938",
                "glory.2791.2791.2",
                "msgfplus_2021_03_22",
                "2791",
                None,
                None,
                "2",
                "40",
                "2.6986221E-12",
                None,
                "Pep_YICDNQDTISSK",
                "4.4458354E-15",
                "3",
                "40",
            ]
        ]
    ).all()


def test_engine_parsers_msgfplus_check_dataframe_integrity_unknown_mod():
    input_file = (
        pytest._test_path / "data" / "BSA1_msgfplus_2021_03_22_unknown_mod.mzid"
    )
    rt_lookup_path = pytest._test_path / "data" / "BSA1_ursgal_lookup.csv"
    db_path = pytest._test_path / "data" / "BSA.fasta"

    parser = MSGFPlus_2021_03_22_Parser(
        input_file,
        params={
            "cpus": 2,
            "rt_pickle_name": rt_lookup_path,
            "database": db_path,
            "enzyme": "(?<=[KR])(?![P])",
            "terminal_cleavage_site_integrity": "any",
            "validation_score_field": {"msgfplus_2021_03_22": "ms-gf:spec_evalue"},
            "bigger_scores_better": {"msgfplus_2021_03_22": False},
            "modifications": [
                {
                    "aa": "M",
                    "type": "opt",
                    "position": "any",
                    "name": "Oxidation",
                },
                {
                    "aa": "C",
                    "type": "opt",
                    "position": "any",
                    "name": "Carbamidomethyl",
                },
                {
                    "aa": "C",
                    "type": "opt",
                    "position": "any",
                    "name": "DTB-IAA",
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
    assert pytest.approx(df["exp_mz"].mean()) == 488.0319
    assert len(df) == 92
    # Currently this excludes two peptides and is reduced
    assert pytest.approx(df["ucalc_mz"].mean()) == 486.5571
    assert (df["raw_data_location"] == "path/for/glory.mzML").all()
    assert df["modifications"].str.contains("Acetyl:0").sum() == 0
    assert df["modifications"].str.contains("Oxidation:").sum() == 0
    # the unknown modification is actually DTB-IAA
    assert df["modifications"].str.contains("DTB-IAA:3").sum() == 2
    assert (
        df[df["sequence"] != "EACFAVEGPK"]["modifications"].str.count(
            "Carbamidomethyl:"
        )
        == df[df["sequence"] != "EACFAVEGPK"]["sequence"].str.count("C")
    ).all()
    assert df["modifications"].str.count(":").sum() == 71
