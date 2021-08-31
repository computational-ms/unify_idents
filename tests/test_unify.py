#!/usr/bin/env python
from pathlib import Path

import pytest
import uparma
from unify_idents.engine_parsers.ident.omssa_parser import OmssaParser
from unify_idents.unify import UnifiedDataFrame, Unify
from unify_idents.engine_parsers.ident.msgfplus_2021_03_22_parser import (
    MSGFPlus_2021_03_22,
)
from unify_idents.engine_parsers.ident.msamanda_parser import MSamandaParser


def test_unify_get_parser_classes():
    rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
    p = Path(__file__).parent / "data" / "test_Creinhardtii_QE_pH11_omssa_2_1_9.csv"
    db_path = Path(__file__).parent / "data" / "BSA.fasta"
    u = Unify(
        p,
        {
            "rt_pickle_name": rt_lookup_path,
            "database": db_path,
            "modifications": [
                "C,fix,any,Carbamidomethyl",
                "M,opt,any,Oxidation",
                "*,opt,Prot-N-term,Acetyl",
            ],
            "omssa_mod_dir": Path(__file__).parent / "data",
        },
    )
    parsers = u._get_parser_classes()
    assert (
        len(parsers) == 8
    )  # currently msamanda, msfragger, msgfplus, omssa, xtandem, flash_lfq and 2 dummies


def test_unify_get_omssa_parser():
    rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
    p = Path(__file__).parent / "data" / "test_Creinhardtii_QE_pH11_omssa_2_1_9.csv"
    db_path = Path(__file__).parent / "data" / "BSA.fasta"
    u = Unify(
        p,
        {
            "rt_pickle_name": rt_lookup_path,
            "database": db_path,
            "modifications": [
                "C,fix,any,Carbamidomethyl",
                "M,opt,any,Oxidation",
                "*,opt,Prot-N-term,Acetyl",
            ],
            "omssa_mod_dir": Path(__file__).parent / "data",
        },
    )
    parser = u._get_parser(p)
    assert isinstance(parser, OmssaParser)


def test_unify_get_msgfplus_parser():
    rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
    p = Path(__file__).parent / "data" / "BSA1_msgfplus_2021_03_22.mzid"
    db_path = Path(__file__).parent / "data" / "BSA.fasta"
    u = Unify(
        p,
        {
            "rt_pickle_name": rt_lookup_path,
            "database": db_path,
            "modifications": [
                "C,fix,any,Carbamidomethyl",
                "M,opt,any,Oxidation",
                "*,opt,Prot-N-term,Acetyl",
            ],
        },
    )
    parser = u._get_parser(p)
    assert isinstance(parser, MSGFPlus_2021_03_22)


def test_engine_parsers_omssa_unified_frame():
    input_file = (
        Path(__file__).parent / "data" / "test_Creinhardtii_QE_pH11_omssa_2_1_9.csv"
    )
    rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
    db_path = Path(__file__).parent / "data" / "BSA.fasta"

    u = Unify(
        input_file,
        {
            "rt_pickle_name": rt_lookup_path,
            "database": db_path,
            "modifications": [
                "C,fix,any,Carbamidomethyl",
                "M,opt,any,Oxidation",
                "*,opt,Prot-N-term,Acetyl",
            ],
            "omssa_mod_dir": Path(__file__).parent / "data",
        },
    )
    df = u.get_dataframe()
    assert isinstance(df, UnifiedDataFrame)


def test_engine_parsers_msfragger_unified_frame():
    input_file = (
        Path(__file__).parent
        / "data"
        / "test_Creinhardtii_QE_pH11_mzml2mgf_0_0_1_msfragger_3.tsv"
    )
    rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
    db_path = Path(__file__).parent / "data" / "test_Creinhardtii_target_decoy.fasta"

    u = Unify(
        input_file,
        params={
            "rt_pickle_name": rt_lookup_path,
            "database": db_path,
            "modifications": [
                "C,fix,any,Carbamidomethyl",
                "M,opt,any,Oxidation",
                "*,opt,Prot-N-term,Acetyl",
            ],
            "omssa_mod_dir": Path(__file__).parent / "data",
            "Raw file location": "test_Creinhardtii_QE_pH11.mzML",
        },
    )
    df = u.get_dataframe()
    assert isinstance(df, UnifiedDataFrame)


def test_engine_parsers_msgf_unified_frame():
    input_file = Path(__file__).parent / "data" / "BSA1_msgfplus_2021_03_22.mzid"
    rt_lookup_path = Path(__file__).parent / "data" / "BSA1_ursgal_lookup.csv.bz2"
    db_path = Path(__file__).parent / "data" / "BSA.fasta"

    u = Unify(
        input_file,
        params={
            "rt_pickle_name": rt_lookup_path,
            "database": db_path,
            "modifications": [
                "C,fix,any,Carbamidomethyl",
                "M,opt,any,Oxidation",
                "*,opt,Prot-N-term,Acetyl",
            ],
        },
    )
    df = u.get_dataframe()
    assert isinstance(df, UnifiedDataFrame)
    assert len(df) == 91
    assert (
        df.df.iloc[0]["Protein ID"]
        == "sp|P02769|ALBU_BOVIN Serum albumin OS=Bos taurus GN=ALB PE=1 SV=4"
    )
    assert df.df.iloc[0]["Sequence Pre AA"] == "K"
    assert df.df.iloc[0]["Sequence Post AA"] == "L"


def test_unify_msfragger_df_masses():
    input_file = Path(__file__).parent / "data" / "msfragger_no_mods.tsv"
    rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
    db_path = Path(__file__).parent / "data" / "test_Creinhardtii_target_decoy.fasta"

    parser = Unify(
        input_file,
        {
            "rt_pickle_name": rt_lookup_path,
            "database": db_path,
            "modifications": [
                "C,fix,any,Carbamidomethyl",
                "M,opt,any,Oxidation",
                "*,opt,Prot-N-term,Acetyl",
            ],
            "Raw file location": "test_Creinhardtii_QE_pH11.mzML",
            "15N": False,
        },
    )
    res = parser.get_dataframe()
    row = res.df.iloc[0]
    assert row["Sequence"] == "ATTALTDDTLDGAGR"
    assert row["Charge"] == "2"
    assert float(row["uCalc m/z"]) == pytest.approx(739.3601)
    assert float(row["uCalc m/z"]) == pytest.approx(739.3601)
    assert float(row["uCalc Mass"]) == pytest.approx(1477.7128)
    assert float(row["Accuracy (ppm)"]) == pytest.approx(-2.182, 0.01)


def test_unify_get_msamanda_parser():
    rt_lookup_path = Path(__file__).parent / "data" / "BSA1_ursgal_lookup.csv.bz2"
    p = Path(__file__).parent / "data" / "BSA_msamanda_2_0_0_17442.csv"
    db_path = Path(__file__).parent / "data" / "BSA.fasta"
    u = Unify(
        p,
        {
            "rt_pickle_name": rt_lookup_path,
            "database": db_path,
            "modifications": [
                "C,fix,any,Carbamidomethyl",
                "M,opt,any,Oxidation",
                "*,opt,Prot-N-term,Acetyl",
            ],
        },
    )
    parser = u._get_parser(p)
    assert isinstance(parser, MSamandaParser)


def test_engine_parsers_msamanda_unified_frame():
    input_file = Path(__file__).parent / "data" / "BSA1_msamanda_2_0_0_17442.csv"
    rt_lookup_path = Path(__file__).parent / "data" / "BSA1_ursgal_lookup.csv.bz2"
    db_path = Path(__file__).parent / "data" / "BSA.fasta"

    u = Unify(
        input_file,
        params={
            "rt_pickle_name": rt_lookup_path,
            "database": db_path,
            "modifications": [
                "C,fix,any,Carbamidomethyl",
                "M,opt,any,Oxidation",
                "*,opt,Prot-N-term,Acetyl",
            ],
        },
    )
    df = u.get_dataframe()
    assert isinstance(df, UnifiedDataFrame)
    assert len(df) == 87
    assert (
        df.df.iloc[0]["Protein ID"]
        == "sp|P02769|ALBU_BOVIN Serum albumin OS=Bos taurus GN=ALB PE=1 SV=4"
    )
    assert df.df.iloc[0]["Sequence Pre AA"] == "K"
    assert df.df.iloc[0]["Sequence Post AA"] == "D"


def test_unify_xtandem_df_masses():
    input_file = (
        Path(__file__).parent / "data" / "test_Creinhardtii_QE_pH11_xtandem_alanine.xml"
    )
    rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
    db_path = Path(__file__).parent / "data" / "test_Creinhardtii_target_decoy.fasta"

    parser = Unify(
        input_file,
        {
            "rt_pickle_name": rt_lookup_path,
            "database": db_path,
            "modifications": [
                "C,fix,any,Carbamidomethyl",
                "M,opt,any,Oxidation",
                "*,opt,Prot-N-term,Acetyl",
            ],
            "Raw file location": "test_Creinhardtii_QE_pH11.mzML",
            "15N": False,
        },
    )
    res = parser.get_dataframe()
    row = res.df[res.df["Sequence"] == "ASDGKYVDEYFAATYVCTDHGRGK"]
    # breakpoint()
    assert row["Sequence"].iloc[0] == "ASDGKYVDEYFAATYVCTDHGRGK"
    assert row["Charge"].iloc[0] == "3"
    assert float(row["uCalc m/z"].iloc[0]) == pytest.approx(
        904.0782, abs=5e-6 * 904.0782
    )
    assert float(row["Calc m/z"].iloc[0]) == pytest.approx(
        904.0782, abs=5e-6 * 904.0782
    )
    assert float(row["uCalc Mass"].iloc[0]) == pytest.approx(
        2710.2202, abs=5e-6 * 2710.2202
    )
    assert row["Modifications"].iloc[0] == "Carbamidomethyl:17"


def test_unify_msgf_df_masses():
    input_file = (
        Path(__file__).parent
        / "data"
        / "test_Creinhardtii_QE_pH11_msgfplus_2021_03_22.mzid"
    )
    rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
    db_path = Path(__file__).parent / "data" / "test_Creinhardtii_target_decoy.fasta"

    parser = Unify(
        input_file,
        {
            "rt_pickle_name": rt_lookup_path,
            "database": db_path,
            "modifications": [
                "C,fix,any,Carbamidomethyl",
                "M,opt,any,Oxidation",
                "*,opt,Prot-N-term,Acetyl",
            ],
            "Raw file location": "test_Creinhardtii_QE_pH11.mzML",
            "15N": False,
        },
    )
    res = parser.get_dataframe()
    row = res.df[res.df["Sequence"] == "ASDGKYVDEYFAATYVCTDHGRGK"]
    assert row["Sequence"].iloc[0] == "ASDGKYVDEYFAATYVCTDHGRGK"
    assert row["Charge"].iloc[0] == "3"
    assert row["Modifications"].iloc[0] == "Carbamidomethyl:17"
    assert float(row["uCalc m/z"].iloc[0]) == pytest.approx(
        904.0782, abs=5e-6 * 904.0782
    )
    assert float(row["Calc m/z"].iloc[0]) == pytest.approx(
        904.0782, abs=5e-6 * 904.0782
    )
    assert float(row["uCalc Mass"].iloc[0]) == pytest.approx(
        2710.2202, abs=5e-6 * 2710.2202
    )


def test_unify_msamanda_df_masses():
    input_file = (
        Path(__file__).parent
        / "data"
        / "test_Creinhardtii_QE_pH11_msamanda_2_0_0_17442.csv"
    )
    rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
    db_path = Path(__file__).parent / "data" / "test_Creinhardtii_target_decoy.fasta"

    parser = Unify(
        input_file,
        {
            "rt_pickle_name": rt_lookup_path,
            "database": db_path,
            "modifications": [
                "C,fix,any,Carbamidomethyl",
                "M,opt,any,Oxidation",
                "*,opt,Prot-N-term,Acetyl",
            ],
            "Raw file location": "test_Creinhardtii_QE_pH11.mzML",
            "15N": False,
        },
    )
    res = parser.get_dataframe()
    row = res.df[res.df["Sequence"] == "ASDGKYVDEYFAATYVCTDHGRGK"]
    assert row["Sequence"].iloc[0] == "ASDGKYVDEYFAATYVCTDHGRGK"
    assert row["Charge"].iloc[0] == "3"
    assert row["Modifications"].iloc[0] == "Carbamidomethyl:17"
    assert float(row["uCalc m/z"].iloc[0]) == pytest.approx(
        904.0782, abs=5e-6 * 904.0782
    )
    assert float(row["Calc m/z"].iloc[0]) == pytest.approx(
        904.0782, abs=5e-6 * 904.0782
    )
    assert float(row["uCalc Mass"].iloc[0]) == pytest.approx(
        2710.2202, abs=5e-6 * 2710.2202
    )
