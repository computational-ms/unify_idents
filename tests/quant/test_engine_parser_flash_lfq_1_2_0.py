#!/usr/bin/env python
from pathlib import Path
from unify_idents.unify import UnifiedDataFrame
from unify_idents.engine_parsers.quant.flash_lfq_1_2_0_parser import FlashLFQ


def test_engine_parsers_flashLFQ_init():
    input_file = (
        Path(__file__).parent.parent / "data" / "flash_lfq_1_2_0_quantified_peaks.tsv"
    )
    rt_lookup_path = (
        Path(__file__).parent.parent / "data" / "BSA1_ursgal_lookup.csv.bz2"
    )
    db_path = Path(__file__).parent / "data" / "BSA.fasta"
    parser = FlashLFQ(input_file, params={"rt_pickle_name": rt_lookup_path})


def test_engine_parsers_flashLFQ_file_matches_parser():
    input_file = (
        Path(__file__).parent.parent / "data" / "flash_lfq_1_2_0_quantified_peaks.tsv"
    )
    rt_lookup_path = (
        Path(__file__).parent.parent / "data" / "BSA1_ursgal_lookup.csv.bz2"
    )
    db_path = Path(__file__).parent.parent / "data" / "BSA.fasta"

    assert FlashLFQ.file_matches_parser(input_file) is True


def test_engine_parsers_flashLFQ_file_not_matches_parser():
    input_file = (
        Path(__file__).parent.parent
        / "data"
        / "test_Creinhardtii_QE_pH11_mzml2mgf_0_0_1_msfragger_3.tsv"
    )
    rt_lookup_path = (
        Path(__file__).parent.parent / "data" / "BSA1_ursgal_lookup.csv.bz2"
    )
    db_path = Path(__file__).parent.parent / "data" / "BSA.fasta"

    assert FlashLFQ.file_matches_parser(input_file) is False


def test_engine_parsers_flashLFQ_iterable():
    input_file = (
        Path(__file__).parent.parent / "data" / "flash_lfq_1_2_0_quantified_peaks.tsv"
    )
    rt_lookup_path = (
        Path(__file__).parent.parent / "data" / "BSA1_ursgal_lookup.csv.bz2"
    )
    db_path = (
        Path(__file__).parent.parent / "data" / "test_Creinhardtii_target_decoy.fasta"
    )

    parser = FlashLFQ(
        input_file,
        params={"rt_pickle_name": rt_lookup_path},
    )
    for row in parser:
        print(row)


def test_engine_parsers_flashLFQ_set_file_name():
    input_file = (
        Path(__file__).parent.parent / "data" / "flash_lfq_1_2_0_quantified_peaks.tsv"
    )
    rt_lookup_path = Path(__file__).parent.parent / "data" / "_ursgal_lookup.csv.bz2"
    db_path = (
        Path(__file__).parent.parent / "data" / "test_Creinhardtii_target_decoy.fasta"
    )

    fname = "/test/path/BSA1.mzML"

    parser = FlashLFQ(
        input_file,
        params={
            "rt_pickle_name": rt_lookup_path,
            "Raw file location": fname,
        },
    )
    for row in parser:
        assert row["trivial_name"] == "YLYEIAR"
        assert row["file_name"] == fname
        assert row["charge"] == "2"
        break


def test_engine_parser_flashLFQ_extract_mods():
    input_file = (
        Path(__file__).parent.parent / "data" / "flash_lfq_1_2_0_quantified_peaks.tsv"
    )
    rt_lookup_path = Path(__file__).parent.parent / "data" / "_ursgal_lookup.csv.bz2"
    db_path = (
        Path(__file__).parent.parent / "data" / "test_Creinhardtii_target_decoy.fasta"
    )

    parser = FlashLFQ(
        input_file,
        params={"rt_pickle_name": rt_lookup_path},
    )
    test_sequence = "ELC[Carbamidomethyl]"
    mods = parser.extract_mods(test_sequence)
    assert mods == "Carbamidomethyl:3"

    test_sequence2 = "ELC[Carbamidomethyl]MMMM[Oxidation]"
    mods = parser.extract_mods(test_sequence2)
    assert mods == "Carbamidomethyl:3;Oxidation:7"

    test_sequence3 = "[Acetyl]ELC[Carbamidomethyl]MMMM[Oxidation]"
    mods = parser.extract_mods(test_sequence3)
    assert mods == "Acetyl:0;Carbamidomethyl:3;Oxidation:7"

    # TODO encode C-terminal mods
    # test_sequence4 = "[Acetyl]ELC[Carbamidomethyl]MMMM[Oxidation][TERMINALMOD]"
    # mods = parser.extract_mods(test_sequence4)
    # assert mods == "Acetyl:0;Carbamidomethyl:3;Oxidation:7;TERMINALMOD:8"
