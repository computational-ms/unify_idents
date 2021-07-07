#!/usr/bin/env python
from pathlib import Path
from unify_idents.unify import UnifiedDataFrame
from unify_idents.engine_parsers.msamanda_parser import MSamandaParser
import uparma

from collections.abc import Iterable


def test_engine_parsers_msamanda_init():
    input_file = (
        Path(__file__).parent / "data" / "BSA1_msamanda_2_0_0_17442.csv"
    )
    rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
    db_path = Path(__file__).parent / "data" / "BSA1.fasta"

    parser = MSamandaParser(
        input_file,
        params={
            "scan_rt_lookup_file": rt_lookup_path,
            "database": db_path,
            "Modifications": [
                "C,fix,any,Carbamidomethyl",
                "M,opt,any,Oxidation",
                "*,opt,Prot-N-term,Acetyl",
            ],
            # "msamanda_mod_dir": Path(__file__).parent / "data",
        },
    )


def test_engine_parsers_msamanda_file_matches_parser():
    input_file = (
            Path(
                __file__).parent / "data" / "BSA1_msamanda_2_0_0_17442.csv"
    )
    rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
    db_path = Path(__file__).parent / "data" / "BSA1.fasta"

    assert MSamandaParser.file_matches_parser(input_file) is True


def test_engine_parsers_msamanda_iterable():
    input_file = (
            Path(
                __file__).parent / "data" / "BSA1_msamanda_2_0_0_17442.csv"
    )
    rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
    db_path = Path(__file__).parent / "data" / "BSA1.fasta"

    parser = MSamandaParser(
        input_file,
        params={
            "scan_rt_lookup_file": rt_lookup_path,
            "database": db_path,
            "Modifications": [
                "C,fix,any,Carbamidomethyl",
                "M,opt,any,Oxidation",
                "*,opt,Prot-N-term,Acetyl",
            ],
            # "Raw file location": "test_Creinhardtii_QE_pH11.mzML",
        },
    )
    # breakpoint()
    assert isinstance(parser, Iterable)


def test_engine_parsers_msamanda_unify_row():
    input_file = (
            Path(
                __file__).parent / "data" / "BSA1_msamanda_2_0_0_17442.csv"
    )
    rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
    db_path = Path(__file__).parent / "data" / "BSA1.fasta"

    parser = MSamandaParser(
        input_file,
        params={
            "scan_rt_lookup_file": rt_lookup_path,
            "database": db_path,
            "Modifications": [
                "C,fix,any,Carbamidomethyl",
                "M,opt,any,Oxidation",
                "*,opt,Prot-N-term,Acetyl",
            ],
            # "omssa_mod_dir": Path(__file__).parent / "data",
        },
    )

    for row in parser:
        print(row)
