from pathlib import Path

from unify_idents.engine_parsers.base_parser import __IdentBaseParser


def test_engine_parsers_IdentBaseParser_init():
    input_file = (
        Path(__file__).parent / "data" / "test_Creinhardtii_QE_pH11_xtandem_alanine.xml"
    )
    rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
    db_path = Path(__file__).parent / "data" / "test_Creinhardtii_target_decoy.fasta"

    parser = __IdentBaseParser(
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
        },
    )


def test_engine_parsers_IdentBaseParser_file_matches_parser_non_existing():
    # should always return False
    __IdentBaseParser.check_parser_compatibility("whatever") is False


def test_engine_parsers_IdentBaseParser_file_matches_parser_existing():
    # should always return False
    __IdentBaseParser.check_parser_compatibility(
        Path(__file__).parent / "data" / "test_Creinhardtii_QE_pH11_xtandem_alanine.xml"
    ) is False
