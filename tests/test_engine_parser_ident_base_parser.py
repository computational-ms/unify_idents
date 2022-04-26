from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from unify_idents.engine_parsers.base_parser import (
    IdentBaseParser,
    get_mass_and_composition,
    merge_and_join_dicts,
)


def test_base_parser_read_rt_lookup_file():
    rt_lookup_path = pytest._test_path / "data" / "BSA1_ursgal_lookup.csv"
    input_file = (
        pytest._test_path / "data" / "test_Creinhardtii_QE_pH11_xtandem_alanine.xml"
    )

    bp = IdentBaseParser(input_file, params={"rt_pickle_name": rt_lookup_path})
    rt_lookup = bp._read_meta_info_lookup_file()
    assert (rt_lookup["rt_unit"] == 1).all()
    assert pytest.approx(
        rt_lookup["precursor_mz"].mean(), 550.8444810049874
    )  # check consistency
    assert pytest.approx(rt_lookup.loc[2450, "precursor_mz"], 618.2697)
    assert pytest.approx(rt_lookup.loc[2450, "rt"], 92067.714)


def test_engine_parsers_IdentBaseParser_init():
    input_file = (
        pytest._test_path / "data" / "test_Creinhardtii_QE_pH11_xtandem_alanine.xml"
    )
    rt_lookup_path = pytest._test_path / "data" / "_ursgal_lookup.csv"
    db_path = pytest._test_path / "data" / "test_Creinhardtii_target_decoy.fasta"

    parser = IdentBaseParser(
        input_file,
        params={
            "cpus": 2,
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


def test_engine_parsers_IdentBaseParser_check_parser_compatibility_non_existing():
    # should always return False
    IdentBaseParser.check_parser_compatibility("whatever") is False


def test_engine_parsers_IdentBaseParser_check_parser_compatibility_existing():
    # should always return False
    IdentBaseParser.check_parser_compatibility(
        pytest._test_path / "data" / "test_Creinhardtii_QE_pH11_xtandem_alanine.xml"
    ) is False


def test_engine_parsers_IdentBase_Parser_sanitize():
    obj = IdentBaseParser(input_file=None, params=None)
    obj.mapping_dict = {"Engine:C": None, "Engine:B": None, "Engine:A": None}
    obj.df = pd.DataFrame(
        np.random.random((5, len(obj.col_order) + 4)),
        columns=obj.col_order.to_list()
        + ["Engine:C", "Engine:B", "Engine:A", "This should not exist"],
    )
    obj.df["raw_data_location"] = "test.mgf"
    obj.df["spectrum_title"] = "spec_title.this_should_go"
    obj.df.loc[3, "raw_data_location"] = None
    obj.df["modifications"] = "ZZ:1;AA:8"
    obj.sanitize()
    assert obj.df.loc[3, "raw_data_location"] == ""
    assert (obj.df["modifications"] == "ZZ:1;AA:8").all()
    assert obj.df.columns.to_list() == obj.col_order.to_list() + [
        "Engine:A",
        "Engine:B",
        "Engine:C",
    ]


def test_add_ranks_increasing_engine_scores_better():
    obj = IdentBaseParser(
        input_file=None,
        params={
            "validation_score_field": {
                "translated_value": {"msfragger_3_0": "msfragger:hyperscore"}
            },
            "bigger_scores_better": {"translated_value": {"msfragger_3_0": True}},
        },
    )
    obj.df = pd.DataFrame(
        np.ones((5, len(obj.col_order) + 1)),
        columns=obj.col_order.to_list() + ["msfragger:hyperscore"],
    )
    obj.df["msfragger:hyperscore"] = [5, 2, 1, 3, 3]
    obj.df["search_engine"] = "msfragger_3_0"
    obj.add_ranks()
    assert obj.df["rank"].to_list() == [1, 4, 5, 2, 2]


def test_add_ranks_decreasing_engine_scores_better():
    obj = IdentBaseParser(
        input_file=None,
        params={
            "validation_score_field": {
                "translated_value": {"msgfplus_2021_03_22": "ms-gf:spec_evalue"}
            },
            "bigger_scores_better": {
                "translated_value": {"msgfplus_2021_03_22": False}
            },
        },
    )
    obj.df = pd.DataFrame(
        np.ones((5, len(obj.col_order) + 1)),
        columns=obj.col_order.to_list() + ["msfragger:hyperscore"],
    )
    obj.df["ms-gf:spec_evalue"] = [5, 2, 1, 3, 3]
    obj.df["search_engine"] = "msgfplus_2021_03_22"
    obj.add_ranks()
    assert obj.df["rank"].to_list() == [5, 2, 1, 3, 3]


def test_add_protein_ids():
    obj = IdentBaseParser(
        input_file=None,
        params={
            "cpus": 2,
            "database": pytest._test_path / "data/test_Creinhardtii_target_decoy.fasta",
            "validation_score_field": {
                "translated_value": {"msfragger_3_0": "msfragger:hyperscore"}
            },
            "bigger_scores_better": {"translated_value": {"msfragger_3_0": True}},
        },
    )
    obj.df = pd.DataFrame(
        np.ones((1, 1)),
        columns=["msfragger:hyperscore"],
    )
    obj.df.loc[0, "sequence"] = "SAVVGTFFR"
    obj.add_protein_ids()

    assert set(obj.df.columns) == {
        "sequence_stop",
        "sequence",
        "protein_id",
        "sequence_post_aa",
        "sequence_pre_aa",
        "sequence_start",
        "msfragger:hyperscore",
    }
    assert obj.df.loc[0, "sequence_start"] == "287"
    assert obj.df.loc[0, "sequence_stop"] == "295"
    assert obj.df.loc[0, "sequence_pre_aa"] == "R"
    assert obj.df.loc[0, "sequence_post_aa"] == "D"
    assert (
        obj.df.loc[0, "protein_id"]
        == "Cre12.g514050.t1.2 pacid=30792611 transcript=Cre12.g514050.t1.2 locus=Cre12.g514050 ID=Cre12.g514050.t1.2.v5.5 annot-version=v5.5"
    )


def test_calc_masses_offsets_and_composition():
    obj = IdentBaseParser(
        input_file=None,
        params={
            "cpus": 2,
            "rt_pickle_name": pytest._test_path / "data/_ursgal_lookup.csv",
            "validation_score_field": {
                "translated_value": {"msfragger_3_0": "msfragger:hyperscore"}
            },
            "bigger_scores_better": {"translated_value": {"msfragger_3_0": True}},
        },
    )
    obj.df = pd.DataFrame(
        np.ones((3, len(obj.col_order) + 1)),
        columns=obj.col_order.to_list() + ["msfragger:hyperscore"],
    )
    obj.df["spectrum_id"] = 3 * [10152]
    obj.get_meta_info()
    obj.df.loc[:, "sequence"] = 3 * ["PEPTCIDE"]
    obj.df.loc[:, "modifications"] = [
        "",
        "Carbamidomethyl:5",
        "Acetyl:0;Carbamidomethyl:5",
    ]
    obj.calc_masses_offsets_and_composition()
    ref_masses = np.array(
        [902.3691, 902.3691 + 57.021464, 902.3691 + 57.021464 + 42.010565]
    )

    assert np.allclose(obj.df["ucalc_mass"], ref_masses, atol=1e-4)
    assert np.allclose(obj.df["ucalc_mz"], ref_masses + obj.PROTON, atol=1e-4)
    assert np.allclose(
        obj.df["accuracy_ppm"], [-159398.095, -209306.942, -242444.593], atol=1e-3
    )


def test_get_exp_rt_and_mz():
    obj = IdentBaseParser(
        input_file=None,
        params={
            "cpus": 2,
            "rt_pickle_name": pytest._test_path / "data/_ursgal_lookup.csv",
        },
    )
    obj.df = pd.DataFrame(
        np.ones((5, len(obj.col_order) + 1)),
        columns=obj.col_order.to_list() + ["msfragger:hyperscore"],
    )
    obj.df["spectrum_id"] = [10152, 10381, 10414, 10581, 11535]
    obj.get_meta_info()
    assert np.allclose(
        obj.df["exp_mz"], [759.379, 439.196, 739.358, 664.286, 496.264], atol=1e-3
    )
    assert np.allclose(
        obj.df["retention_time_seconds"],
        [114735.00, 116583.52, 116805.30, 118033.46, 124967.99],
        atol=1e-2,
    )


def test_create_mod_dicts():
    obj = IdentBaseParser(
        input_file=None,
        params={
            "cpus": 2,
            "modifications": [
                {
                    "aa": "*",
                    "type": "opt",
                    "position": "Prot-N-term",
                    "name": "Acetyl",
                },
                {
                    "aa": "C",
                    "type": "fix",
                    "position": "any",
                    "name": "Carbamidomethyl",
                },
            ],
        },
    )
    mod_dict = obj._create_mod_dicts()
    reference_dict = {
        "Carbamidomethyl": {"mass": 57.021464, "aa": {"C", "any"}, "position": {"any"}},
        "Acetyl": {
            "mass": 42.010565,
            "aa": {"Prot-N-term", "*"},
            "position": {"Prot-N-term"},
        },
    }
    assert mod_dict == reference_dict


def test_calc_mz():
    obj = IdentBaseParser(input_file=None, params=None)
    masses_in_weird_types = pd.Series(["10", 20, 30.0])
    charges_in_weird_types = pd.Series(["2", 2, 2.0])
    mz = obj._calc_mz(masses_in_weird_types, charges_in_weird_types).to_list()

    assert mz[0] == pytest.approx(6.0, abs=1e-2)
    assert mz[1] == pytest.approx(11.0, abs=1e-2)
    assert mz[2] == pytest.approx(16.0, abs=1e-2)


def test_get_mass_and_composition():
    seq = "PEPTCIDE"
    mods = ""
    # Without modifications
    m, comp = get_mass_and_composition(seq=seq, mods=mods)
    assert m == pytest.approx(902.3691)
    assert comp == "C(37)H(58)N(8)O(16)S(1)"

    # With modifications
    mods = "Acetyl:0;Carbamidomethyl:5"
    m, comp = get_mass_and_composition(seq=seq, mods=mods)
    assert m == pytest.approx(902.3691 + 57.021464 + 42.010565)
    assert comp == "C(41)H(63)N(9)O(18)S(1)"


def test_merge_and_join_dicts():
    dict_a = {"a": "part_a", "b": "part_a"}
    dict_b = {"a": "part_b", "b": "part_b"}
    out_dict = merge_and_join_dicts([dict_a, dict_b], ";")
    assert out_dict == {"a": "part_a;part_b", "b": "part_a;part_b"}


def test_assert_only_iupac_and_missing_aas():
    obj = IdentBaseParser(
        input_file=None,
        params={
            "cpus": 2,
            "rt_pickle_name": pytest._test_path / "data/_ursgal_lookup.csv",
        },
    )
    obj.df = pd.DataFrame(
        np.ones((4, len(obj.col_order) + 1)),
        columns=obj.col_order.to_list() + ["msfragger:hyperscore"],
    )
    obj.df["sequence"] = ["PEPTIDE", "[3jd2]", "ACDXU", "MREPEPTIDE"]
    obj.assert_only_iupac_and_missing_aas()

    assert all(obj.df.index == [0, 3])


def test_add_decoy_identity():
    obj = IdentBaseParser(
        input_file=None,
        params={
            "cpus": 2,
            "rt_pickle_name": pytest._test_path / "data/_ursgal_lookup.csv",
        },
    )
    obj.df = pd.DataFrame(
        np.ones((4, len(obj.col_order) + 1)),
        columns=obj.col_order.to_list() + ["msfragger:hyperscore"],
    )
    obj.df["protein_id"] = ["NOTADECOY", "PEPTIDE", "decoy_PEPTIDE", "decoy_ASDF"]

    obj.add_decoy_identity()

    assert all(obj.df["is_decoy"] == [False, False, True, True])


def test_add_decoy_identity_non_default_prefix():
    obj = IdentBaseParser(
        input_file=None,
        params={
            "cpus": 2,
            "rt_pickle_name": pytest._test_path / "data/_ursgal_lookup.csv",
            "decoy_tag": "non_default_tag_",
        },
    )
    obj.df = pd.DataFrame(
        np.ones((4, len(obj.col_order) + 1)),
        columns=obj.col_order.to_list() + ["msfragger:hyperscore"],
    )
    obj.df["protein_id"] = [
        "NOTADECOY",
        "PEPTIDE",
        "decoy_PEPTIDE",
        "non_default_tag_ASDF",
    ]

    obj.add_decoy_identity()

    assert all(obj.df["is_decoy"] == [False, False, False, True])


def test_check_enzyme_specificity_trypsin_all():
    obj = IdentBaseParser(
        input_file=None,
        params={
            "cpus": 2,
            "rt_pickle_name": pytest._test_path / "data/_ursgal_lookup.csv",
            "enzyme": {
                "original_value": "trypsin",
                "translated_value": "(?<=[KR])(?![P])",
            },
            "terminal_cleavage_site_integrity": {"translated_value": "all"},
        },
    )
    obj.df = pd.DataFrame(
        np.ones((4, len(obj.col_order) + 1)),
        columns=obj.col_order.to_list() + ["msfragger:hyperscore"],
    )
    obj.df["sequence"] = ["PEPRTIDEK", "EPTIDEK", "EPTIDEK", "EPRPTRIRDEK"]
    obj.df["sequence_pre_aa"] = ["K", "K", "A", "K"]
    obj.df["sequence_post_aa"] = ["A", "P", "A<|>P", "A<|>V"]
    obj.check_enzyme_specificity()

    assert all(obj.df["enzn"] == [False, True, False, True])
    assert all(obj.df["enzc"] == [True, False, False, True])
    assert all(obj.df["missed_cleavages"] == [1, 0, 0, 2])


def test_check_enzyme_specificity_trypsin_any():
    obj = IdentBaseParser(
        input_file=None,
        params={
            "cpus": 2,
            "rt_pickle_name": pytest._test_path / "data/_ursgal_lookup.csv",
            "enzyme": {
                "original_value": "trypsin",
                "translated_value": "(?<=[KR])(?![P])",
            },
            "terminal_cleavage_site_integrity": {"translated_value": "any"},
        },
    )
    obj.df = pd.DataFrame(
        np.ones((4, len(obj.col_order) + 1)),
        columns=obj.col_order.to_list() + ["msfragger:hyperscore"],
    )
    obj.df["sequence"] = ["PEPRTIDEK", "EPTIDEK", "EPTIDEK", "EPTIDEK"]
    obj.df["sequence_pre_aa"] = ["K", "K", "A", "K"]
    obj.df["sequence_post_aa"] = ["A", "P", "A<|>P", "A<|>V"]
    obj.check_enzyme_specificity()

    assert all(obj.df["enzn"] == [False, True, False, True])
    assert all(obj.df["enzc"] == [True, False, True, True])
    assert all(obj.df["missed_cleavages"] == [1, 0, 0, 0])


def test_check_enzyme_specificity_nonspecific():
    obj = IdentBaseParser(
        input_file=None,
        params={
            "cpus": 2,
            "rt_pickle_name": pytest._test_path / "data/_ursgal_lookup.csv",
            "enzyme": {"original_value": "nonspecific", "translated_value": ".^"},
            "terminal_cleavage_site_integrity": {"translated_value": "all"},
        },
    )
    obj.df = pd.DataFrame(
        np.ones((4, len(obj.col_order) + 1)),
        columns=obj.col_order.to_list() + ["msfragger:hyperscore"],
    )
    obj.df["sequence"] = ["ASDF", "EPTIDEK", "EPTIDEK", "EPTIDEK"]
    obj.df["sequence_pre_aa"] = ["L<|>E", "A", "R", "P"]
    obj.df["sequence_post_aa"] = ["F", "L<|>E", "A<|>P", "A<|>V"]
    obj.check_enzyme_specificity()

    assert all(obj.df["enzn"] == [True, True, True, True])
    assert all(obj.df["enzc"] == [True, True, True, True])
    assert all(obj.df["missed_cleavages"] == [0, 0, 0, 0])
