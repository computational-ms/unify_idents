#!/usr/bin/env python
from pathlib import Path
from unify_idents.unify import UnifiedDataFrame
from unify_idents.engine_parsers.ident.msgfplus_2021_03_22_parser import (
    MSGFPlus_2021_03_22,
)
import uparma
from collections.abc import Iterable
from unify_idents import UnifiedRow


def test_engine_parsers_msgfplus_init():
    input_file = Path(__file__).parent.parent / "data" / "BSA1_msgfplus_2021_03_22.mzid"
    rt_lookup_path = Path(__file__).parent.parent / "data" / "BSA1_ursgal_lookup.csv.bz2"
    db_path = (
        Path(__file__).parent.parent / "data" / "test_Creinhardtii_target_decoy.fasta"
    )

    parser = MSGFPlus_2021_03_22(
        input_file,
        params={
            "rt_pickle_name": rt_lookup_path,
            "database": db_path,
            "modifications": [
                {
                    "aa": "M",
                    "type": "opt",
                    "position": "any",
                    "name": "Oxidation",
                },
                {
                    "aa": "C",
                    "type": "fix",
                    "position": "any",
                    "name": "Carbamidomethyl",
                },
                {
                    "aa": "*",
                    "type": "opt",
                    "position": "Prot-N-term",
                    "name": "Acetyl",
                },
            ],
            "omssa_mod_dir": Path(__file__).parent.parent / "data",
        },
    )


def test_engine_parsers_msgfplus_file_matches_parser():
    msgf_parser_class = MSGFPlus_2021_03_22
    input_file = Path(__file__).parent.parent / "data" / "BSA1_msgfplus_2021_03_22.mzid"
    assert msgf_parser_class.file_matches_parser(input_file) is True


def test_engine_parsers_msgfplus_file_matches_parser_fail_with_omsa_file():
    msgf_parser_class = MSGFPlus_2021_03_22
    input_file = (
        Path(__file__).parent.parent
        / "data"
        / "test_Creinhardtii_QE_pH11_omssa_2_1_9.csv"
    )
    assert msgf_parser_class.file_matches_parser(input_file) is False


def test_engine_parsers_msgfplus_is_iterable():
    input_file = Path(__file__).parent.parent / "data" / "BSA1_msgfplus_2021_03_22.mzid"
    rt_lookup_path = Path(__file__).parent.parent / "data" / "BSA1_ursgal_lookup.csv.bz2"
    db_path = (
        Path(__file__).parent.parent / "data" / "test_Creinhardtii_target_decoy.fasta"
    )

    parser = MSGFPlus_2021_03_22(
        input_file,
        params={
            "rt_pickle_name": rt_lookup_path,
            "database": db_path,
            "modifications": [
                {
                    "aa": "M",
                    "type": "opt",
                    "position": "any",
                    "name": "Oxidation",
                },
                {
                    "aa": "C",
                    "type": "fix",
                    "position": "any",
                    "name": "Carbamidomethyl",
                },
                {
                    "aa": "*",
                    "type": "opt",
                    "position": "Prot-N-term",
                    "name": "Acetyl",
                },
            ],
            "omssa_mod_dir": Path(__file__).parent.parent / "data",
        },
    )

    assert isinstance(parser, Iterable)


def test_engine_parsers_msgfplus_iter_items():
    input_file = Path(__file__).parent.parent / "data" / "BSA1_msgfplus_2021_03_22.mzid"
    rt_lookup_path = Path(__file__).parent.parent / "data" / "BSA1_ursgal_lookup.csv.bz2"
    db_path = (
        Path(__file__).parent.parent / "data" / "test_Creinhardtii_target_decoy.fasta"
    )

    parser = MSGFPlus_2021_03_22(
        input_file,
        params={
            "rt_pickle_name": rt_lookup_path,
            "database": db_path,
            "modifications": [
                {
                    "aa": "M",
                    "type": "opt",
                    "position": "any",
                    "name": "Oxidation",
                },
                {
                    "aa": "C",
                    "type": "fix",
                    "position": "any",
                    "name": "Carbamidomethyl",
                },
                {
                    "aa": "*",
                    "type": "opt",
                    "position": "Prot-N-term",
                    "name": "Acetyl",
                },
            ],
            "omssa_mod_dir": Path(__file__).parent.parent / "data",
        },
    )

    for i, line in enumerate(parser):
        assert isinstance(line, UnifiedRow)


def test_engine_parsers_msgfplus_get_peptide_lookup():
    input_file = Path(__file__).parent.parent / "data" / "BSA1_msgfplus_2021_03_22.mzid"
    rt_lookup_path = Path(__file__).parent.parent / "data" / "BSA1_ursgal_lookup.csv.bz2"
    db_path = (
        Path(__file__).parent.parent / "data" / "test_Creinhardtii_target_decoy.fasta"
    )

    parser = MSGFPlus_2021_03_22(
        input_file,
        params={
            "rt_pickle_name": rt_lookup_path,
            "database": db_path,
            "modifications": [
                {
                    "aa": "M",
                    "type": "opt",
                    "position": "any",
                    "name": "Oxidation",
                },
                {
                    "aa": "C",
                    "type": "fix",
                    "position": "any",
                    "name": "Carbamidomethyl",
                },
                {
                    "aa": "*",
                    "type": "opt",
                    "position": "Prot-N-term",
                    "name": "Acetyl",
                },
            ],
            "omssa_mod_dir": Path(__file__).parent.parent / "data",
        },
    )
    assert len(parser.peptide_lookup) == 24
    assert "Pep_YICDNQDTISSK" in parser.peptide_lookup.keys()
    assert parser.peptide_lookup["Pep_YICDNQDTISSK"]["Sequence"] == "YICDNQDTISSK"
    assert parser.peptide_lookup["Pep_YICDNQDTISSK"]["Modifications"][0]["pos"] == "3"
    assert (
        parser.peptide_lookup["Pep_YICDNQDTISSK"]["Modifications"][0]["mass"]
        == "57.021464"
    )
    assert (
        parser.peptide_lookup["Pep_YICDNQDTISSK"]["Modifications"][0]["name"]
        == "Carbamidomethyl"
    )


def test_engine_parsers_msgfplus_internal_next():
    input_file = Path(__file__).parent.parent / "data" / "BSA1_msgfplus_2021_03_22.mzid"
    rt_lookup_path = Path(__file__).parent.parent / "data" / "BSA1_ursgal_lookup.csv.bz2"
    db_path = (
        Path(__file__).parent.parent / "data" / "test_Creinhardtii_target_decoy.fasta"
    )

    parser = MSGFPlus_2021_03_22(
        input_file,
        params={
            "rt_pickle_name": rt_lookup_path,
            "database": db_path,
            "modifications": [
                {
                    "aa": "M",
                    "type": "opt",
                    "position": "any",
                    "name": "Oxidation",
                },
                {
                    "aa": "C",
                    "type": "fix",
                    "position": "any",
                    "name": "Carbamidomethyl",
                },
                {
                    "aa": "*",
                    "type": "opt",
                    "position": "Prot-N-term",
                    "name": "Acetyl",
                },
            ],
            "omssa_mod_dir": Path(__file__).parent.parent / "data",
        },
    )
    for row in parser._next():
        assert isinstance(row, dict)
        assert (
            row["Raw data location"]
            == "/Users/cellzome/Dev/Gits/Ursgal/ursgal_master/example_data/BSA1.mgf"
        )
        assert row["Peptide"] == "YICDNQDTISSK"
        assert row["Modifications"] == "Carbamidomethyl:3"
        assert row["MS-GF:RawScore"] == "40"
        assert row["MS-GF:NumMatchedMainIons"] == "3"
        # search engine not set here, just stuff directly from the mzid
        break


def test_engine_parsers_msgfplus_next():
    input_file = Path(__file__).parent.parent / "data" / "BSA1_msgfplus_2021_03_22.mzid"
    rt_lookup_path = Path(__file__).parent.parent / "data" / "BSA1_ursgal_lookup.csv.bz2"
    db_path = (
        Path(__file__).parent.parent / "data" / "test_Creinhardtii_target_decoy.fasta"
    )

    parser = MSGFPlus_2021_03_22(
        input_file,
        params={
            "rt_pickle_name": rt_lookup_path,
            "database": db_path,
            "modifications": [
                {
                    "aa": "M",
                    "type": "opt",
                    "position": "any",
                    "name": "Oxidation",
                },
                {
                    "aa": "C",
                    "type": "fix",
                    "position": "any",
                    "name": "Carbamidomethyl",
                },
                {
                    "aa": "*",
                    "type": "opt",
                    "position": "Prot-N-term",
                    "name": "Acetyl",
                },
            ],
            "omssa_mod_dir": Path(__file__).parent.parent / "data",
        },
    )
    row = next(parser)
    assert isinstance(row, UnifiedRow)
    assert row["Sequence"] == "YICDNQDTISSK"
    assert row["Modifications"] == "Carbamidomethyl:3"
    assert row["MS-GF:RawScore"] == "40"
    assert row["MS-GF:NumMatchedMainIons"] == "3"
    assert row["Search Engine"] == "msgfplus_2021_03_22"


def test_engine_parsers_msgfplus_multiple_psms():
    input_file = Path(__file__).parent.parent / "data" / "BSA1_msgfplus_2021_03_22.mzid"
    rt_lookup_path = Path(__file__).parent.parent / "data" / "BSA1_ursgal_lookup.csv.bz2"
    db_path = (
        Path(__file__).parent.parent / "data" / "test_Creinhardtii_target_decoy.fasta"
    )

    parser = MSGFPlus_2021_03_22(
        input_file,
        params={
            "rt_pickle_name": rt_lookup_path,
            "database": db_path,
            "modifications": [
                {
                    "aa": "M",
                    "type": "opt",
                    "position": "any",
                    "name": "Oxidation",
                },
                {
                    "aa": "C",
                    "type": "fix",
                    "position": "any",
                    "name": "Carbamidomethyl",
                },
                {
                    "aa": "*",
                    "type": "opt",
                    "position": "Prot-N-term",
                    "name": "Acetyl",
                },
            ],
            "omssa_mod_dir": Path(__file__).parent.parent / "data",
        },
    )
    for i, row in enumerate(parser):
        pass
    assert i == 91
