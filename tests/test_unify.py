#!/usr/bin/env python
from pathlib import Path
import uparma

from unify_idents.unify import Unify
from unify_idents.engine_parsers.omssa_parser import OmssaParser
from unify_idents.unify import UnifiedDataFrame


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
    assert len(parsers) == 3  # currently we have dummy and omssa


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
