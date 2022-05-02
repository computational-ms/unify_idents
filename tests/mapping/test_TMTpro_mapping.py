import pandas as pd

COMET_TEST_FILE = ""
MASCOT_TEST_FILE = ""
MSAMANDA_TEST_FILE = ""
MSFRAGGER_TEST_FILE = ""
MSGFPLUS_TEST_FILE = ""
OMSSA_TEST_FILE = ""
XTANDEM_TEST_FILE = ""


def test_comet_mapping():
    df = pd.read_csv(COMET_TEST_FILE)
    assert df["modifications"].str.contains("TMTpro").all() is True
    assert df["modifications"].str.contains("Carbamidomethyl").all() is True


def test_mascot_mapping():
    df = pd.read_csv(MASCOT_TEST_FILE)
    assert df["modifications"].str.contains("TMTpro").all() is True
    assert df["modifications"].str.contains("Carbamidomethyl").all() is True


def test_msamanda_mapping():
    df = pd.read_csv(MSAMANDA_TEST_FILE)
    assert df["modifications"].str.contains("TMTpro").all() is True
    assert df["modifications"].str.contains("Carbamidomethyl").all() is True


def test_msfragger_mapping():
    df = pd.read_csv(MSFRAGGER_TEST_FILE)
    assert df["modifications"].str.contains("TMTpro").all() is True
    assert df["modifications"].str.contains("Carbamidomethyl").all() is True


def test_msgfplus_mapping():
    df = pd.read_csv(MSGFPLUS_TEST_FILE)
    assert df["modifications"].str.contains("TMTpro").all() is True
    assert df["modifications"].str.contains("Carbamidomethyl").all() is True


def test_omssa_mapping():
    df = pd.read_csv(OMSSA_TEST_FILE)
    assert df["modifications"].str.contains("TMTpro").all() is True
    assert df["modifications"].str.contains("Carbamidomethyl").all() is True


def test_xtandem_mapping():
    df = pd.read_csv(XTANDEM_TEST_FILE)
    assert df["modifications"].str.contains("TMTpro").all() is True
    assert df["modifications"].str.contains("Carbamidomethyl").all() is True
