#!/usr/bin/env python
from pathlib import Path
from unify_idents.unify import UnifiedDataFrame
from unify_idents.engine_parsers.omssa_parser import OmssaParser
import uparma


def test_engine_parsers_omssa_init():
    input_file = (
        Path(__file__).parent / "data" / "BSA1_mzml2mgf_0_0_1_omssa_2_1_9.csv_tmp"
    )
    rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
    db_path = Path(__file__).parent / "data" / "BSA1.fasta"

    parser = OmssaParser(
        input_file,
        params={
            "scan_rt_lookup_file": rt_lookup_path,
            "database": db_path,
            "Modifications": [
                "C,fix,any,Carbamidomethyl",
                "M,opt,any,Oxidation",
            ],
        },
    )


def test_engine_parsers_omssa_unify_row():
    input_file = (
        Path(__file__).parent / "data" / "BSA1_mzml2mgf_0_0_1_omssa_2_1_9.csv_tmp"
    )
    rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
    db_path = Path(__file__).parent / "data" / "BSA1.fasta"

    parser = OmssaParser(
        input_file,
        params={
            "scan_rt_lookup_file": rt_lookup_path,
            "database": db_path,
            "Modifications": [
                "C,fix,any,Carbamidomethyl",
                "M,opt,any,Oxidation",
            ],
        },
    )
    for row in parser:
        print(row)


# def test_engine_parsers_omssa_unified_frame():
#     input_file = (
#         Path(__file__).parent / "data" / "BSA1_mzml2mgf_0_0_1_omssa_2_1_9.csv_tmp"
#     )
#     rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
#     db_path = Path(__file__).parent / "data" / "BSA1.fasta"

#     parser = OmssaParser(
#         input_file,
#         params={
#             "scan_rt_lookup_file": rt_lookup_path,
#             "database": db_path,
#             "Modifications": [
#                 "C,fix,any,Carbamidomethyl",
#                 "M,opt,any,Oxidation",
#             ],
#         },
#     )
#     df = parser.get_dataframe()
#     assert isinstance(df, UnifiedDataFrame)
