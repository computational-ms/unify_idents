#!/usr/bin/env python
from pathlib import Path

import pytest
import uparma
from unify_idents.engine_parsers.omssa_parser import OmssaParser
from unify_idents.unify import UnifiedDataFrame, Unify
from unify_idents.engine_parsers.msgfplus_2021_03_22_parser import MSGFPlus_2021_03_22


def test_unify_get_parser_classes():
    rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
    p = Path(__file__).parent / "data" / "test_Creinhardtii_QE_pH11_omssa_2_1_9.csv"
    db_path = Path(__file__).parent / "data" / "test_Creinhardtii_target_decoy.fasta"
    u = Unify(
        p,
        {
            "scan_rt_lookup_file": rt_lookup_path,
            "database": db_path,
            "Modifications": [
                "C,fix,any,Carbamidomethyl",
                "M,opt,any,Oxidation",
                "*,opt,Prot-N-term,Acetyl",
            ],
            "omssa_mod_dir": Path(__file__).parent / "data",
        },
    )
    parsers = u._get_parser_classes()
    assert len(parsers) == 4  # currently we have dummy and omssa


def test_unify_get_omssa_parser():
    rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
    p = Path(__file__).parent / "data" / "test_Creinhardtii_QE_pH11_omssa_2_1_9.csv"
    db_path = Path(__file__).parent / "data" / "test_Creinhardtii_target_decoy.fasta"
    u = Unify(
        p,
        {
            "scan_rt_lookup_file": rt_lookup_path,
            "database": db_path,
            "Modifications": [
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
    db_path = Path(__file__).parent / "data" / "test_Creinhardtii_target_decoy.fasta"
    u = Unify(
        p,
        {
            "scan_rt_lookup_file": rt_lookup_path,
            "database": db_path,
            "Modifications": [
                "C,fix,any,Carbamidomethyl",
                "M,opt,any,Oxidation",
                "*,opt,Prot-N-term,Acetyl",
            ],
            "omssa_mod_dir": Path(__file__).parent / "data",
        },
    )
    parser = u._get_parser(p)
    assert isinstance(parser, MSGFPlus_2021_03_22)


def test_engine_parsers_omssa_unified_frame():
    input_file = (
        Path(__file__).parent / "data" / "test_Creinhardtii_QE_pH11_omssa_2_1_9.csv"
    )
    rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
    db_path = Path(__file__).parent / "data" / "test_Creinhardtii_target_decoy.fasta"

    u = Unify(
        input_file,
        {
            "scan_rt_lookup_file": rt_lookup_path,
            "database": db_path,
            "Modifications": [
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
            "scan_rt_lookup_file": rt_lookup_path,
            "database": db_path,
            "Modifications": [
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
    rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
    db_path = Path(__file__).parent / "data" / "test_Creinhardtii_target_decoy.fasta"

    u = Unify(
        input_file,
        params={
            "scan_rt_lookup_file": rt_lookup_path,
            "database": db_path,
            "Modifications": [
                "C,fix,any,Carbamidomethyl",
                "M,opt,any,Oxidation",
                "*,opt,Prot-N-term,Acetyl",
            ],
            "omssa_mod_dir": Path(__file__).parent / "data",
        },
    )
    df = u.get_dataframe()
    assert isinstance(df, UnifiedDataFrame)


def test_unify_msfragger_df_masses():
    input_file = Path(__file__).parent / "data" / "msfragger_no_mods.tsv"
    rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
    db_path = Path(__file__).parent / "data" / "test_Creinhardtii_target_decoy.fasta"

    parser = Unify(
        input_file,
        {
            "scan_rt_lookup_file": rt_lookup_path,
            "database": db_path,
            "Modifications": [
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
    assert float(row["uCalc Mass"]) == pytest.approx(1476.7056)
    assert float(row["Accuracy (ppm)"]) == pytest.approx(-2.182, 0.01)
