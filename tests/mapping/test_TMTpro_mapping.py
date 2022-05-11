import pytest

import unify_idents

COMET_TEST_FILE = pytest._test_path / "data" / "mapping_data" / "trunc_comet.mzid"
# MASCOT_TEST_FILE = ""
# MSAMANDA_TEST_FILE = ""
MSFRAGGER_TEST_FILE = pytest._test_path / "data" / "mapping_data" / "trunc_msfragger.tsv"
MSGFPLUS_TEST_FILE = pytest._test_path / "data" / "mapping_data" / "trunc_msgfplus.mzid"
OMSSA_TEST_FILE = pytest._test_path / "data" / "mapping_data" / "trunc_omssa.csv"
XTANDEM_TEST_FILE = pytest._test_path / "data" / "mapping_data" / "trunc_xtandem.xml"

PARAMS = {
    "database": pytest._test_path / "data" / "mapping_data" / "trunc_fasta.protein.faa",
    "rt_pickle_name": pytest._test_path / "data" / "mapping_data" / "trunc_meta.spectra_meta.csv",
    "enzyme": {
        "original_value": "trypsin",
        "translated_value": "(?<=[KR])(?![P])",
    },
    "terminal_cleavage_site_integrity": {"translated_value": "any"},
    "validation_score_field": {
        "translated_value": {
            "msfragger_3_0": "msfragger:hyperscore",
            "comet_2020_01_4": "comet:e_value",
            "msgfplus_2021_03_22": "ms-gf:spec_evalue",
            "xtandem_alanine": "x!tandem:hyperscore",
            "omssa_2_1_9": "omssa:pvalue",
        }
    },
    "bigger_scores_better": {
        "translated_value": {
            "msfragger_3_0": True,
            "comet_2020_01_4": False,
            "msgfplus_2021_03_22": False,
            "xtandem_alanine": True,
            "omssa_2_1_9": True,
        }
    },
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
        {
            "aa": "K",
            "type": "fix",
            "position": "any",
            "name": "TMTpro",
        },
        {
            "aa": "*",
            "type": "opt",
            "position": "N-term",
            "name": "TMTpro",
        },
    ],
}


def test_comet_mapping():
    df = unify_idents.Unify(input_file=COMET_TEST_FILE, params=PARAMS).get_dataframe()
    mod_str = df.loc[0, "modifications"]
    assert mod_str == "TMTpro:0;Carbamidomethyl:8;Oxidation:12;Carbamidomethyl:19"


# def test_mascot_mapping():
#     df = pd.read_csv(MASCOT_TEST_FILE)
#     assert df["modifications"].str.contains("TMTpro").all() is True
#     assert df["modifications"].str.contains("Carbamidomethyl").all() is True
#
#
# def test_msamanda_mapping():
#     df = pd.read_csv(MSAMANDA_TEST_FILE)
#     assert df["modifications"].str.contains("TMTpro").all() is True
#     assert df["modifications"].str.contains("Carbamidomethyl").all() is True


def test_msfragger_mapping():
    df = unify_idents.Unify(
        input_file=MSFRAGGER_TEST_FILE, params=PARAMS
    ).get_dataframe()
    mod_str = df.loc[0, "modifications"]
    assert mod_str == "TMTpro:0;Carbamidomethyl:6"


def test_msgfplus_mapping():
    df = unify_idents.Unify(
        input_file=MSGFPLUS_TEST_FILE, params=PARAMS
    ).get_dataframe()
    mod_str = df.loc[0, "modifications"]
    assert mod_str == "TMTpro:0;Oxidation:7;Carbamidomethyl:14;TMTpro:26;TMTpro:27"


def test_omssa_mapping():
    df = unify_idents.Unify(
        input_file=OMSSA_TEST_FILE, params=PARAMS
    ).get_dataframe()
    mod_str = df.loc[0, "modifications"]
    assert mod_str == "TMTpro:0;Carbamidomethyl:6"


def test_xtandem_mapping():
    df = unify_idents.Unify(input_file=XTANDEM_TEST_FILE, params=PARAMS).get_dataframe()
    mod_str = df.loc[0, "modifications"]
    assert mod_str == "TMTpro:0;Oxidation:6;Carbamidomethyl:7"
