from unify_idents.engine_parsers.msfragger3_parser import MSFragger3Parser
from pathlib import Path


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


def test_engine_parsers_msfragger_file_matches_parser():
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
            "scan_rt_lookup_file": rt_lookup_path,
            "database": db_path,
            "Modifications": [
                "C,fix,any,Carbamidomethyl",
                "M,opt,any,Oxidation",
                "*,opt,Prot-N-term,Acetyl",
            ],
            "Raw file location": "test_Creinhardtii_QE_pH11.mzML",
            # "omssa_mod_dir": Path(__file__).parent / "data",
        },
    )
    assert parser.file_matches_parser() is True


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
            "scan_rt_lookup_file": rt_lookup_path,
            "database": db_path,
            "Modifications": [
                "C,fix,any,Carbamidomethyl",
                "M,opt,any,Oxidation",
                "*,opt,Prot-N-term,Acetyl",
            ],
            "Raw file location": "test_Creinhardtii_QE_pH11.mzML",
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
            "scan_rt_lookup_file": rt_lookup_path,
            "database": db_path,
            "Modifications": [
                "C,fix,any,Carbamidomethyl",
                "M,opt,any,Oxidation",
                "*,opt,Prot-N-term,Acetyl",
            ],
            "Raw file location": "test_Creinhardtii_QE_pH11.mzML",
        },
    )
    for row in parser:
        print(row)
        assert row["Protein ID"] == ""
        assert row["Sequence"] == ""
        assert row["uCalc m/z"] == 0
        break
