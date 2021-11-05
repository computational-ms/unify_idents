#!/usr/bin/env python
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
