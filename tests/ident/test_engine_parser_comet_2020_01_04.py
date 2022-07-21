#!/usr/bin/env python

import pandas as pd
import pytest
from lxml import etree

from unify_idents.engine_parsers.ident.comet_2020_01_4_parser import (
    Comet_2020_01_4_Parser,
    _get_single_spec_df,
)


def test_engine_parsers_comet_init():
    input_file = pytest._test_path / "data" / "BSA1_comet_2020_01_4.mzid"

    parser = Comet_2020_01_4_Parser(
        input_file,
        params={
            "cpus": 2,
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
    input_file = (
        pytest._test_path / "data" / "test_Creinhardtii_QE_pH11_omssa_2_1_9.csv"
    )
    assert msgf_parser_class.check_parser_compatibility(input_file) is False


def test_engine_parsers_comet_check_dataframe_integrity():
    input_file = pytest._test_path / "data" / "BSA1_comet_2020_01_4.mzid"
    rt_lookup_path = pytest._test_path / "data" / "BSA1_ursgal_lookup.csv"
    db_path = pytest._test_path / "data" / "BSA.fasta"

    parser = Comet_2020_01_4_Parser(
        input_file,
        params={
            "cpus": 2,
            "enzyme": "(?<=[KR])(?![P])",
            "terminal_cleavage_site_integrity": "any",
            "validation_score_field": {"comet_2020_01_4": "comet:e_value"},
            "bigger_scores_better": {"comet_2020_01_4": False},
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
    assert pytest.approx(df["ucalc_mz"].mean()) == 457.85944
    assert pytest.approx(df["exp_mz"].mean()) == 457.87625

    assert df["modifications"].str.contains("Acetyl:0").sum() == 5
    assert df["modifications"].str.contains("Oxidation:").sum() == 0
    assert (
        df["modifications"].str.count("Carbamidomethyl:")
        == df["sequence"].str.count("C")
    ).all()
    assert df["modifications"].str.count(":").sum() == 38
    assert (df["raw_data_location"] == "path/for/glory.mzML").all()


def test_get_single_spec_df():
    input_file = pytest._test_path / "data" / "BSA1_comet_2020_01_4.mzid"
    element = (
        etree.parse(input_file)
        .getroot()
        .find(".//{*}SpectrumIdentificationList/{*}SpectrumIdentificationResult")
    )
    ref_dict = {
        "exp_mz": None,
        "calc_mz": None,
        "spectrum_title": None,
        "search_engine": "comet_2020_01_4",
        "spectrum_id": None,
        "modifications": None,
        "retention_time_seconds": None,
        "charge": None,
        "comet:score": None,
        "comet:deltacn": None,
        "comet:xcorr": None,
        "comet:evalue": None,
        "sequence": None,
        "comet:spec_evalue": None,
        "comet:num_matched_ions": None,
        "comet:num_unmatched_ions": None,
    }
    mapping_dict = {
        "chargeState": "charge",
        "Comet:spscore": "comet:score",
        "Comet:deltacn": "comet:deltacn",
        "Comet:xcorr": "comet:xcorr",
        "Comet:expectation value": "comet:evalue",
        "peptide_ref": "sequence",
        "experimentalMassToCharge": "exp_mz",
        "calculatedMassToCharge": "calc_mz",
        "SpecEValue": "comet:spec_evalue",
        "number of matched peaks": "comet:num_matched_ions",
        "number of unmatched peaks": "comet:num_unmatched_ions",
    }
    _get_single_spec_df.reference_dict = ref_dict
    _get_single_spec_df.mapping_dict = mapping_dict
    result = _get_single_spec_df(etree.tostring(element))

    assert isinstance(result, pd.DataFrame)
    assert (
        result.values
        == [
            [
                "358.174682",
                "358.174575",
                None,
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
