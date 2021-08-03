#!/usr/bin/env python
from pathlib import Path
from unify_idents.unify import UnifiedDataFrame

from unify_idents.engine_parsers.ident.omssa_parser import OmssaParser
import uparma


def test_engine_parsers_omssa_init():
    input_file = (
        Path(__file__).parent.parent
        / "data"
        / "test_Creinhardtii_QE_pH11_omssa_2_1_9.csv"
    )
    rt_lookup_path = Path(__file__).parent.parent / "data" / "_ursgal_lookup.csv.bz2"
    db_path = (
        Path(__file__).parent.parent / "data" / "test_Creinhardtii_target_decoy.fasta"
    )

    parser = OmssaParser(
        input_file,
        params={
            "rt_pickle_name": rt_lookup_path,
            "database": db_path,
            "modifications": [
                "C,fix,any,Carbamidomethyl",
                "M,opt,any,Oxidation",
                "*,opt,Prot-N-term,Acetyl",
            ],
            "omssa_mod_dir": Path(__file__).parent.parent / "data",
        },
    )


def test_engine_parsers_omssa_file_matches_parser():
    input_file = (
        Path(__file__).parent.parent
        / "data"
        / "test_Creinhardtii_QE_pH11_omssa_2_1_9.csv"
    )
    rt_lookup_path = Path(__file__).parent.parent / "data" / "_ursgal_lookup.csv.bz2"
    db_path = (
        Path(__file__).parent.parent / "data" / "test_Creinhardtii_target_decoy.fasta"
    )

    assert OmssaParser.file_matches_parser(input_file) is True


def test_engine_parsers_omssa_unify_row():
    input_file = (
        Path(__file__).parent.parent
        / "data"
        / "test_Creinhardtii_QE_pH11_omssa_2_1_9.csv"
    )
    rt_lookup_path = Path(__file__).parent.parent / "data" / "_ursgal_lookup.csv.bz2"
    db_path = (
        Path(__file__).parent.parent / "data" / "test_Creinhardtii_target_decoy.fasta"
    )

    parser = OmssaParser(
        input_file,
        params={
            "rt_pickle_name": rt_lookup_path,
            "database": db_path,
            "modifications": [
                "C,fix,any,Carbamidomethyl",
                "M,opt,any,Oxidation",
                "*,opt,Prot-N-term,Acetyl",
            ],
            "omssa_mod_dir": Path(__file__).parent.parent / "data",
        },
    )
    for row in parser:
        print(row)


# def test_engine_parsers_omssa_unified_frame():
#     input_file = (
#         Path(__file__).parent.parent / "data" / "BSA1_mzml2mgf_0_0_1_omssa_2_1_9.csv_tmp"
#     )
#     rt_lookup_path = Path(__file__).parent.parent / "data" / "_ursgal_lookup.csv.bz2"
#     db_path = Path(__file__).parent.parent / "data" / "test_Creinhardtii_target_decoy.fasta"

#     parser = OmssaParser(
#         input_file,
#         params={
#             "rt_pickle_name": rt_lookup_path,
#             "database": db_path,
#             "modifications": [
#                 "C,fix,any,Carbamidomethyl",
#                 "M,opt,any,Oxidation",
#             ],
#         },
#     )
#     df = parser.get_dataframe()
#     assert isinstance(df, UnifiedDataFrame)
