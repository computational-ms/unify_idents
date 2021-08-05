#!/usr/bin/env python
from pathlib import Path
from unify_idents.unify import UnifiedDataFrame, UnifiedRow
from unify_idents.engine_parsers.omssa_parser import OmssaParser
import uparma
import pytest


def test_engine_parsers_omssa_init():
    input_file = (
        Path(__file__).parent / "data" / "test_Creinhardtii_QE_pH11_omssa_2_1_9.csv"
    )
    rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
    db_path = Path(__file__).parent / "data" / "test_Creinhardtii_target_decoy.fasta"

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
            "omssa_mod_dir": Path(__file__).parent / "data",
        },
    )


def test_engine_parsers_omssa_file_matches_parser():
    input_file = (
        Path(__file__).parent / "data" / "test_Creinhardtii_QE_pH11_omssa_2_1_9.csv"
    )
    rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
    db_path = Path(__file__).parent / "data" / "test_Creinhardtii_target_decoy.fasta"

    assert OmssaParser.file_matches_parser(input_file) is True


def test_engine_parsers_omssa_unify_row():
    input_file = (
        Path(__file__).parent / "data" / "test_Creinhardtii_QE_pH11_omssa_2_1_9.csv"
    )
    rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
    db_path = Path(__file__).parent / "data" / "test_Creinhardtii_target_decoy.fasta"

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
            "omssa_mod_dir": Path(__file__).parent / "data",
        },
    )
    for row in parser:
        print(row)


def test_engine_parsers_omssa_is_iterable():
    input_file = (
        Path(__file__).parent / "data" / "test_Creinhardtii_QE_pH11_omssa_2_1_9.csv"
    )
    rt_lookup_path = Path(__file__).parent / "data" / "BSA1_ursgal_lookup.csv.bz2"
    db_path = Path(__file__).parent / "data" / "test_Creinhardtii_target_decoy.fasta"

    parser = MSGFPlus_2021_03_22(
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
        },
    )
    assert isinstance(parser, Iterable)



def test_engine_parsers_omssa_next():
    input_file = (
        Path(__file__).parent / "data" / "test_Creinhardtii_QE_pH11_omssa_2_1_9.csv"
    )
    rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
    db_path = Path(__file__).parent / "data" / "test_Creinhardtii_target_decoy.fasta"
    parser = OmssaParser(
    row = next(parser)
    assert isinstance(row, UnifiedRow)
    print(row.data.keys())
    assert row["Sequence"] == "ALAMEWGPFPRLMVVACNDAINVCRK"
    assert row["Modifications"] == "Oxidation:4;Carbamidomethyl:17;Carbamidomethyl:24"
    assert (
        row["Raw data location"]
        == "/Users/cellzome/Dev/Gits/Ursgal/ursgal_master/example_data/BSA1.mgf"
    )
    assert row["Charge"] == "4"
    assert float(row["OMSSA:pvalue"]) == pytest.approx(0.000166970504409832)
    assert float(row["uCalc m/z"]) == 0
    assert float(row["uCalc mass"]) == 0
    assert float(row["Calc mass"]) == pytest.approx(3033.491)
    assert float(row["Calc m/z"]) == pytest.approx(759.38002646677)
    # assert row["MS-GF:RawScore"] == "40"
    # assert row["MS-GF:NumMatchedMainIons"] == "3"
    # assert row["Search Engine"] == "MSGFPlus_2021_03_22"

