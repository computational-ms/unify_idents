from pathlib import Path

import pytest

from unify_idents.engine_parsers.ident.msfragger_3_parser import MSFragger_3_Parser


def test_engine_parsers_msfragger_init():
    input_file = (
        pytest._test_path / "data" / "test_Creinhardtii_QE_pH11_msfragger_3.tsv"
    )
    parser = MSFragger_3_Parser(
        input_file,
        params={
            "cpus": 2,
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
            # "omssa_mod_dir": pytest._test_path / "data",
        },
    )


def test_engine_parsers_msfragger_check_parser_compatibility():
    input_file = (
        pytest._test_path / "data" / "test_Creinhardtii_QE_pH11_msfragger_3.tsv"
    )
    assert MSFragger_3_Parser.check_parser_compatibility(input_file) is True


def test_engine_parsers_msfragger_check_dataframe_integrity():
    input_file = pytest._test_path / "data" / "BSA1_msfragger_3.tsv"
    rt_lookup_path = pytest._test_path / "data" / "BSA1_ursgal_lookup.csv.bz2"
    db_path = pytest._test_path / "data" / "BSA.fasta"

    parser = MSFragger_3_Parser(
        input_file,
        params={
            "cpus": 2,
            "rt_pickle_name": rt_lookup_path,
            "database": db_path,
            "Raw data location": "path/for/glory.mzML",
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
            "Raw data location": "path/for/glory.mzML",
            "15N": False,
        },
    )
    df = parser.unify()
    assert len(df) == 3417
    assert pytest.approx(df["uCalc m/z"].mean()) == 477.8585
    assert pytest.approx(df["Exp m/z"].mean()) == 478.12137

    assert df["Modifications"].str.contains("Acetyl:0").sum() == 2
    assert df["Modifications"].str.contains("Oxidation:").sum() == 221
    assert (
        df["Modifications"].str.count("Carbamidomethyl:")
        == df["Sequence"].str.count("C")
    ).all()
    assert df["Modifications"].str.count(":").sum() == 2242
    assert (df["Raw data location"] == "path/for/glory.mzML").all()


def test_map_mod_translation():
    input_file = (
        pytest._test_path / "data" / "test_Creinhardtii_QE_pH11_msfragger_3.tsv"
    )

    parser = MSFragger_3_Parser(
        input_file,
        params={
            "cpus": 2,
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
            "Raw data location": "path/for/glory.mzML",
            "15N": False,
        },
    )
    map_dict = {"15.994915": ["Oxidation"], "57.021465": ["Carbamidomethyl"]}
    converted = parser._map_mod_translation(
        row=["3M(15.994915)", "15M(15.994915)", "18C(57.021465)"], map_dict=map_dict
    )
    assert converted == "Oxidation:3;Oxidation:15;Carbamidomethyl:18;"
