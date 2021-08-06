from unify_idents.engine_parsers.msfragger3_parser import MSFragger3Parser
from pathlib import Path

from collections import Iterable


def test_engine_parsers_msfragger_init():
    input_file = (
        Path(__file__).parent
        / "data"
        / "test_Creinhardtii_QE_pH11_mzml2mgf_0_0_1_msfragger_3.tsv"
    )
    rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
    db_path = Path(__file__).parent / "data" / "test_Creinhardtii_target_decoy.fasta"

    parser = MSFragger3Parser(
        input_file,
        params={
            "rt_pickle_name": rt_lookup_path,
            "database": db_path,
            "modifications": [
                "C,fix,any,Carbamidomethyl",
                "M,opt,any,Oxidation",
                "*,opt,Prot-N-term,Acetyl",
            ],
            # "omssa_mod_dir": Path(__file__).parent / "data",
        },
    )


def test_engine_parsers_msfragger_file_matches_parser():
    input_file = (
        Path(__file__).parent
        / "data"
        / "test_Creinhardtii_QE_pH11_mzml2mgf_0_0_1_msfragger_3.tsv"
    )
    rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
    db_path = Path(__file__).parent / "data" / "test_Creinhardtii_target_decoy.fasta"
    assert MSFragger3Parser.file_matches_parser(input_file) is True


def test_engine_parsers_msfragger_iterable():
    input_file = (
        Path(__file__).parent
        / "data"
        / "test_Creinhardtii_QE_pH11_mzml2mgf_0_0_1_msfragger_3.tsv"
    )
    rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
    db_path = Path(__file__).parent / "data" / "test_Creinhardtii_target_decoy.fasta"

    parser = MSFragger3Parser(
        input_file,
        params={
            "rt_pickle_name": rt_lookup_path,
            "database": db_path,
            "modifications": [
                "C,fix,any,Carbamidomethyl",
                "M,opt,any,Oxidation",
                "*,opt,Prot-N-term,Acetyl",
            ],
            "Raw data location": "/Users/cellzome/Dev/Gits/Ursgal/ursgal_master/example_data/test_Creinhardtii_QE_pH11.mzML",
            "15N": False,
        },
    )
    for row in parser:
        print(row)


def test_engine_parsers_msfragger_unify_row():
    input_file = (
        Path(__file__).parent
        / "data"
        / "test_Creinhardtii_QE_pH11_mzml2mgf_0_0_1_msfragger_3.tsv"
    )
    rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
    db_path = Path(__file__).parent / "data" / "test_Creinhardtii_target_decoy.fasta"

    parser = MSFragger3Parser(
        input_file,
        params={
            "rt_pickle_name": rt_lookup_path,
            "database": db_path,
            "modifications": [
                "C,fix,any,Carbamidomethyl",
                "M,opt,any,Oxidation",
                "*,opt,Prot-N-term,Acetyl",
            ],
            "Raw data location": "/Users/cellzome/Dev/Gits/Ursgal/ursgal_master/example_data/test_Creinhardtii_QE_pH11.mzML",
            "15N": False,
        },
    )
    for row in parser:
        assert row["Sequence"] == "ATTALTDDTLDGAGR"
        assert row["Search Engine"] == "msfragger_3_0"
        break


def test_engine_parsers_msfragger_merge_mods():
    input_file = Path(__file__).parent / "data" / "msfragger_merged_mods.tsv"
    rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
    db_path = Path(__file__).parent / "data" / "test_Creinhardtii_target_decoy.fasta"

    parser = MSFragger3Parser(
        input_file,
        params={
            "rt_pickle_name": rt_lookup_path,
            "database": db_path,
            "modifications": [
                "C,fix,any,Carbamidomethyl",
                "M,opt,any,Oxidation",
                "*,opt,Prot-N-term,Acetyl",
            ],
            "Raw data location": "/Users/cellzome/Dev/Gits/Ursgal/ursgal_master/example_data/test_Creinhardtii_QE_pH11.mzML",
            "15N": False,
        },
    )
    for line in parser:
        assert line["Modifications"] == "Acetyl:0;Carbamidomethyl:1"
        assert line["Sequence"] == "CGFSTVGSGFGSR"
        assert (
            line["Raw data location"]
            == "/Users/cellzome/Dev/Gits/Ursgal/ursgal_master/example_data/test_Creinhardtii_QE_pH11.mzML"
        )
        # assert float(line["uCalc m/z"]) == 631.2851


def test_engine_parsers_msfragger_single_mods():
    input_file = Path(__file__).parent / "data" / "msfragger_single_mod.tsv"
    rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
    db_path = Path(__file__).parent / "data" / "test_Creinhardtii_target_decoy.fasta"

    parser = MSFragger3Parser(
        input_file,
        params={
            "rt_pickle_name": rt_lookup_path,
            "database": db_path,
            "modifications": [
                "C,fix,any,Carbamidomethyl",
                "M,opt,any,Oxidation",
                "*,opt,Prot-N-term,Acetyl",
            ],
            "Raw data location": "/Users/cellzome/Dev/Gits/Ursgal/ursgal_master/example_data/test_Creinhardtii_QE_pH11.mzML",
            "15N": False,
        },
    )
    for line in parser:
        assert "Carbamidomethyl:1" == line["Modifications"]


def test_engine_parsers_msfragger_single_mods():
    input_file = Path(__file__).parent / "data" / "msfragger_no_mods.tsv"
    rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
    db_path = Path(__file__).parent / "data" / "test_Creinhardtii_target_decoy.fasta"

    parser = MSFragger3Parser(
        input_file,
        params={
            "rt_pickle_name": rt_lookup_path,
            "database": db_path,
            "modifications": [
                "C,fix,any,Carbamidomethyl",
                "M,opt,any,Oxidation",
                "*,opt,Prot-N-term,Acetyl",
            ],
            "Raw data location": "/Users/cellzome/Dev/Gits/Ursgal/ursgal_master/example_data/test_Creinhardtii_QE_pH11.mzML",
            "15N": False,
        },
    )
    for line in parser:
        assert "" == line["Modifications"]
