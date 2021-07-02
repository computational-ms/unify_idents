#!/usr/bin/env python
from pathlib import Path
from unify_idents.unify import UnifiedDataFrame
from unify_idents.engine_parsers.msgfplus_2021_03_22_parser import MSGFPlus_2021_03_22
import uparma
from collections import Iterable
from unify_idents import UnifiedRow


def test_engine_parsers_msgfplus_init():
    input_file = Path(__file__).parent / "data" / "BSA1_msgfplus_2021_03_22.mzid"
    rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
    db_path = Path(__file__).parent / "data" / "test_Creinhardtii_target_decoy.fasta"

    parser = MSGFPlus_2021_03_22(
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


def test_engine_parsers_msgfplus_file_matches_parser():
    msgf_parser_class = MSGFPlus_2021_03_22
    input_file = Path(__file__).parent / "data" / "BSA1_msgfplus_2021_03_22.mzid"
    assert msgf_parser_class.file_matches_parser(input_file) is True


def test_engine_parsers_msgfplus_file_matches_parser_fail_with_omssa_file():
    msgf_parser_class = MSGFPlus_2021_03_22
    input_file = (
        Path(__file__).parent / "data" / "test_Creinhardtii_QE_pH11_omssa_2_1_9.csv"
    )
    assert msgf_parser_class.file_matches_parser(input_file) is False


def test_engine_parsers_msgfplus_is_iterable():
    input_file = Path(__file__).parent / "data" / "BSA1_msgfplus_2021_03_22.mzid"
    rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
    db_path = Path(__file__).parent / "data" / "test_Creinhardtii_target_decoy.fasta"

    parser = MSGFPlus_2021_03_22(
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

    assert isinstance(parser, Iterable)


def test_engine_parsers_msgfplus_iter_items():
    input_file = Path(__file__).parent / "data" / "BSA1_msgfplus_2021_03_22.mzid"
    rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
    db_path = Path(__file__).parent / "data" / "test_Creinhardtii_target_decoy.fasta"

    parser = MSGFPlus_2021_03_22(
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

    for line in parser:
        assert isinstance(line, UnifiedRow)
