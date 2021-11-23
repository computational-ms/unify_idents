from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from unify_idents.engine_parsers.base_parser import (
    IdentBaseParser,
    get_mass_and_composition,
    merge_and_join_dicts,
)


def test_engine_parsers_IdentBaseParser_init():
    input_file = (
        Path(__file__).parent / "data" / "test_Creinhardtii_QE_pH11_xtandem_alanine.xml"
    )
    rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
    db_path = Path(__file__).parent / "data" / "test_Creinhardtii_target_decoy.fasta"

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
        Path(__file__).parent / "data" / "test_Creinhardtii_QE_pH11_xtandem_alanine.xml"
    ) is False


def test_engine_parsers_IdentBase_Parser_sanitize():
    obj = IdentBaseParser(input_file=None, params=None)
    obj.mapping_dict = {"Engine:C": None, "Engine:B": None, "Engine:A": None}
    obj.df = pd.DataFrame(
        np.ones((5, len(obj.col_order) + 4)),
        columns=obj.col_order.to_list()
        + ["Engine:C", "Engine:B", "Engine:A", "This should not exist"],
    )
    obj.df["Raw data location"] = "test.mgf"
    obj.df["Spectrum Title"] = "spec_title.this_should_go"
    obj.df.loc[3, "Raw data location"] = None
    obj.df["Modifications"] = "ZZ:1;AA:8"
    obj.df.drop(columns="Mass Difference", inplace=True)
    obj.sanitize()
    assert not obj.df["Raw data location"].str.contains(".mgf").any()
    assert obj.df.loc[3, "Raw data location"] == "spec_title"
    assert (obj.df["Modifications"] == "ZZ:1;AA:8").all()
    assert obj.df.columns.to_list() == obj.col_order.to_list() + [
        "Engine:A",
        "Engine:B",
        "Engine:C",
    ]


def test_add_ranks_increasing_engine_scores_better():
    obj = IdentBaseParser(input_file=None, params=None)
    obj.df = pd.DataFrame(
        np.ones((5, len(obj.col_order) + 1)),
        columns=obj.col_order.to_list() + ["MSFragger:Hyperscore"],
    )
    obj.df["MSFragger:Hyperscore"] = [5, 2, 1, 3, 3]
    obj.df["Search Engine"] = "msfragger_3_0"
    obj.add_ranks()
    assert obj.df["Rank"].to_list() == [1, 4, 5, 2, 2]


def test_add_ranks_decreasing_engine_scores_better():
    obj = IdentBaseParser(input_file=None, params=None)
    obj.df = pd.DataFrame(
        np.ones((5, len(obj.col_order) + 1)),
        columns=obj.col_order.to_list() + ["MSFragger:Hyperscore"],
    )
    obj.df["MS-GF:SpecEValue"] = [5, 2, 1, 3, 3]
    obj.df["Search Engine"] = "msgfplus_2021_03_22"
    obj.add_ranks()
    assert obj.df["Rank"].to_list() == [5, 2, 1, 3, 3]


def test_add_protein_ids():
    obj = IdentBaseParser(
        input_file=None,
        params={
            "cpus": 2,
            "database": Path(__file__).parent
            / "data/test_Creinhardtii_target_decoy.fasta",
        },
    )
    obj.df = pd.DataFrame(
        np.ones((1, 1)),
        columns=["MSFragger:Hyperscore"],
    )
    obj.df.loc[0, "Sequence"] = "SAVVGTFFR"
    obj.add_protein_ids()

    assert set(obj.df.columns) == {
        "Sequence Stop",
        "Sequence",
        "Protein ID",
        "Sequence Post AA",
        "Sequence Pre AA",
        "Sequence Start",
        "MSFragger:Hyperscore",
    }
    assert obj.df.loc[0, "Sequence Start"] == "287"
    assert obj.df.loc[0, "Sequence Stop"] == "295"
    assert obj.df.loc[0, "Sequence Pre AA"] == "R"
    assert obj.df.loc[0, "Sequence Post AA"] == "D"
    assert (
        obj.df.loc[0, "Protein ID"]
        == "Cre12.g514050.t1.2 pacid=30792611 transcript=Cre12.g514050.t1.2 locus=Cre12.g514050 ID=Cre12.g514050.t1.2.v5.5 annot-version=v5.5"
    )


def test_calc_masses_offsets_and_composition():
    obj = IdentBaseParser(
        input_file=None,
        params={
            "cpus": 2,
            "rt_pickle_name": Path(__file__).parent / "data/_ursgal_lookup.csv.bz2",
        },
    )
    obj.df = pd.DataFrame(
        np.ones((3, len(obj.col_order) + 1)),
        columns=obj.col_order.to_list() + ["MSFragger:Hyperscore"],
    )
    obj.df["Spectrum ID"] = 3 * [10152]
    obj.get_exp_rt_and_mz()
    obj.df.loc[:, "Sequence"] = 3 * ["PEPTCIDE"]
    obj.df.loc[:, "Modifications"] = [
        "",
        "Carbamidomethyl:5",
        "Acetyl:0;Carbamidomethyl:5",
    ]
    obj.calc_masses_offsets_and_composition()
    ref_masses = np.array(
        [902.3691, 902.3691 + 57.021464, 902.3691 + 57.021464 + 42.010565]
    )

    assert np.allclose(obj.df["uCalc Mass"], ref_masses, atol=1e-4)
    assert np.allclose(obj.df["uCalc m/z"], ref_masses + obj.PROTON, atol=1e-4)
    assert np.allclose(
        obj.df["Accuracy (ppm)"], [-159398.095, -209306.942, -242444.593], atol=1e-3
    )


def test_get_exp_rt_and_mz():
    obj = IdentBaseParser(
        input_file=None,
        params={
            "cpus": 2,
            "rt_pickle_name": Path(__file__).parent / "data/_ursgal_lookup.csv.bz2",
        },
    )
    obj.df = pd.DataFrame(
        np.ones((5, len(obj.col_order) + 1)),
        columns=obj.col_order.to_list() + ["MSFragger:Hyperscore"],
    )
    obj.df["Spectrum ID"] = [10152, 10381, 10414, 10581, 11535]
    obj.get_exp_rt_and_mz()
    assert np.allclose(
        obj.df["Exp m/z"], [759.379, 439.196, 739.358, 664.286, 496.264], atol=1e-3
    )
    assert np.allclose(
        obj.df["Retention Time (s)"],
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
            "rt_pickle_name": Path(__file__).parent / "data/_ursgal_lookup.csv.bz2",
        },
    )
    obj.df = pd.DataFrame(
        np.ones((4, len(obj.col_order) + 1)),
        columns=obj.col_order.to_list() + ["MSFragger:Hyperscore"],
    )
    obj.df["Sequence"] = ["PEPTIDE", "[3jd2]", "ACDXU", "MREPEPTIDE"]
    obj.assert_only_iupac_and_missing_aas()

    assert all(obj.df.index == [0, 3])


def test_add_decoy_identity():
    obj = IdentBaseParser(
        input_file=None,
        params={
            "cpus": 2,
            "rt_pickle_name": Path(__file__).parent / "data/_ursgal_lookup.csv.bz2",
        },
    )
    obj.df = pd.DataFrame(
        np.ones((4, len(obj.col_order) + 1)),
        columns=obj.col_order.to_list() + ["MSFragger:Hyperscore"],
    )
    obj.df["Protein ID"] = ["NOTADECOY", "PEPTIDE", "decoy_PEPTIDE", "decoy_ASDF"]

    obj.add_decoy_identity()

    assert all(obj.df["Is decoy"] == [False, False, True, True])


def test_add_decoy_identity_non_default_prefix():
    obj = IdentBaseParser(
        input_file=None,
        params={
            "cpus": 2,
            "rt_pickle_name": Path(__file__).parent / "data/_ursgal_lookup.csv.bz2",
            "decoy_tag": "non_default_tag_",
        },
    )
    obj.df = pd.DataFrame(
        np.ones((4, len(obj.col_order) + 1)),
        columns=obj.col_order.to_list() + ["MSFragger:Hyperscore"],
    )
    obj.df["Protein ID"] = [
        "NOTADECOY",
        "PEPTIDE",
        "decoy_PEPTIDE",
        "non_default_tag_ASDF",
    ]

    obj.add_decoy_identity()

    assert all(obj.df["Is decoy"] == [False, False, False, True])


def test_check_enzyme_specificity_trypsin_all():
    obj = IdentBaseParser(
        input_file=None,
        params={
            "cpus": 2,
            "rt_pickle_name": Path(__file__).parent / "data/_ursgal_lookup.csv.bz2",
            "enzyme": "trypsin",
            "terminal_cleavage_site_integrity": "all",
        },
    )
    obj.df = pd.DataFrame(
        np.ones((4, len(obj.col_order) + 1)),
        columns=obj.col_order.to_list() + ["MSFragger:Hyperscore"],
    )
    obj.df["Sequence"] = ["PEPRTIDEK", "EPTIDEK", "EPTIDEK", "EPRPTRIRDEK"]
    obj.df["Sequence Pre AA"] = ["K", "K", "A", "K"]
    obj.df["Sequence Post AA"] = ["A", "P", "A<|>P", "A<|>V"]
    obj.check_enzyme_specificity()

    assert all(obj.df["enzN"] == [False, True, False, True])
    assert all(obj.df["enzC"] == [True, False, False, True])
    assert all(obj.df["Missed Cleavages"] == [1, 0, 0, 2])


def test_check_enzyme_specificity_trypsin_any():
    obj = IdentBaseParser(
        input_file=None,
        params={
            "cpus": 2,
            "rt_pickle_name": Path(__file__).parent / "data/_ursgal_lookup.csv.bz2",
            "enzyme": "trypsin",
            "terminal_cleavage_site_integrity": "any",
        },
    )
    obj.df = pd.DataFrame(
        np.ones((4, len(obj.col_order) + 1)),
        columns=obj.col_order.to_list() + ["MSFragger:Hyperscore"],
    )
    obj.df["Sequence"] = ["PEPRTIDEK", "EPTIDEK", "EPTIDEK", "EPTIDEK"]
    obj.df["Sequence Pre AA"] = ["K", "K", "A", "K"]
    obj.df["Sequence Post AA"] = ["A", "P", "A<|>P", "A<|>V"]
    obj.check_enzyme_specificity()

    assert all(obj.df["enzN"] == [False, True, False, True])
    assert all(obj.df["enzC"] == [True, False, True, True])
    assert all(obj.df["Missed Cleavages"] == [1, 0, 0, 0])


def test_check_enzyme_specificity_nonspecific():
    obj = IdentBaseParser(
        input_file=None,
        params={
            "cpus": 2,
            "rt_pickle_name": Path(__file__).parent / "data/_ursgal_lookup.csv.bz2",
            "enzyme": "nonspecific",
            "terminal_cleavage_site_integrity": "all",
        },
    )
    obj.df = pd.DataFrame(
        np.ones((4, len(obj.col_order) + 1)),
        columns=obj.col_order.to_list() + ["MSFragger:Hyperscore"],
    )
    obj.df["Sequence"] = ["ASDF", "EPTIDEK", "EPTIDEK", "EPTIDEK"]
    obj.df["Sequence Pre AA"] = ["L<|>E", "A", "R", "P"]
    obj.df["Sequence Post AA"] = ["F", "L<|>E", "A<|>P", "A<|>V"]
    obj.check_enzyme_specificity()

    assert all(obj.df["enzN"] == [True, True, True, True])
    assert all(obj.df["enzC"] == [True, True, True, True])
    assert all(obj.df["Missed Cleavages"] == [0, 0, 0, 0])
