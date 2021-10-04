#!/usr/bin/env python
from unify_idents.engine_parsers.base_parser import __QuantBaseParser
from pathlib import Path


def test_engine_parsers_QuantBaseParser_init():
    input_file = (
        Path(__file__).parent / "data" / "test_Creinhardtii_QE_pH11_xtandem_alanine.xml"
    )
    rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
    db_path = Path(__file__).parent / "data" / "test_Creinhardtii_target_decoy.fasta"

    parser = __QuantBaseParser(
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


def test_engine_parsers_QuantBaseParser_file_matches_parser_non_existing():
    # should always return False
    __QuantBaseParser.file_matches_parser("whatever") is False


def test_engine_parsers_QuantBaseParser_file_matches_parser_existing():
    # should always return False
    __QuantBaseParser.file_matches_parser(
        Path(__file__).parent / "data" / "test_Creinhardtii_QE_pH11_xtandem_alanine.xml"
    ) is False


def test_engine_parsers_QuantBaseParser_check_required_headers_empty_row():
    # should always return False
    rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
    db_path = Path(__file__).parent / "data" / "test_Creinhardtii_target_decoy.fasta"
    qbp = __QuantBaseParser(
        Path(__file__).parent
        / "data"
        / "test_Creinhardtii_QE_pH11_xtandem_alanine.xml",
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
    qbp.check_required_headers({}) is False


def test_engine_parsers_QuantBaseParser_check_required_headers_wrong_col_names():
    # should always return False
    db_path = Path(__file__).parent / "data" / "test_Creinhardtii_target_decoy.fasta"
    rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
    qbp = __QuantBaseParser(
        Path(__file__).parent
        / "data"
        / "test_Creinhardtii_QE_pH11_xtandem_alanine.xml",
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
    qbp.check_required_headers({"dies": "das", "ana": "nas"}) is False


def test_engine_parsers_QuantBaseParser_check_required_headers_all_present():
    # should always return False
    db_path = Path(__file__).parent / "data" / "test_Creinhardtii_target_decoy.fasta"
    rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
    qbp = __QuantBaseParser(
        Path(__file__).parent
        / "data"
        / "test_Creinhardtii_QE_pH11_xtandem_alanine.xml",
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
    row = {key: "" for key in qbp.required_headers}
    qbp.check_required_headers(row) is False


def test_engine_parsers_QuantBaseParser_check_required_headers_all_present():
    # should always return False
    db_path = Path(__file__).parent / "data" / "test_Creinhardtii_target_decoy.fasta"
    rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
    qbp = __QuantBaseParser(
        Path(__file__).parent
        / "data"
        / "test_Creinhardtii_QE_pH11_xtandem_alanine.xml",
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
    row = {key: "" for key in qbp.required_headers}
    row.update({"haette": "haette", "Fahrrad": "kette"})
    qbp.check_required_headers(row) is False
