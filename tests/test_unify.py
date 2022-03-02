#!/usr/bin/env python
from pathlib import Path
import pytest

from unify_idents.engine_parsers.ident.comet_2020_01_4_parser import (
    Comet_2020_01_4_Parser,
)
from unify_idents.engine_parsers.ident.msamanda_2_parser import MSAmanda_2_Parser
from unify_idents.engine_parsers.ident.msfragger_3_parser import MSFragger_3_Parser
from unify_idents.engine_parsers.ident.msgfplus_2021_03_22_parser import (
    MSGFPlus_2021_03_22_Parser,
)
from unify_idents.engine_parsers.ident.omssa_2_1_9_parser import Omssa_Parser
from unify_idents.engine_parsers.ident.xtandem_alanine import XTandemAlanine_Parser
from unify_idents.engine_parsers.ident.mascot_2_6_2_parser import Mascot_2_6_2_Parser
from unify_idents.unify import Unify
import unify_idents


def test_unify_get_parser_classes():
    rt_lookup_path = pytest._test_path / "data" / "_ursgal_lookup.csv.bz2"
    p = pytest._test_path / "data" / "test_Creinhardtii_QE_pH11_omssa_2_1_9.csv"
    db_path = pytest._test_path / "data" / "BSA.fasta"
    u = Unify(
        p,
        {
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
            "omssa_mod_dir": pytest._test_path / "data",
        },
    )
    # Get files and subtract __init__.py
    ident_files = (
        len(
            list(
                (Path(unify_idents.__path__[0]) / "engine_parsers" / "ident").glob(
                    "*.py"
                )
            )
        )
        - 1
    )
    quant_files = (
        len(
            list(
                (Path(unify_idents.__path__[0]) / "engine_parsers" / "quant").glob(
                    "*.py"
                )
            )
        )
        - 1
    )

    assert len(u._parser_classes) == ident_files + quant_files


def test_unify_get_omssa_parser():
    rt_lookup_path = pytest._test_path / "data" / "_ursgal_lookup.csv.bz2"
    p = pytest._test_path / "data" / "test_Creinhardtii_QE_pH11_omssa_2_1_9.csv"
    db_path = pytest._test_path / "data" / "BSA.fasta"
    u = Unify(
        p,
        {
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
            "omssa_mod_dir": pytest._test_path / "data",
        },
    )
    assert isinstance(u.parser, Omssa_Parser)


def test_unify_get_msgfplus_parser():
    rt_lookup_path = pytest._test_path / "data" / "_ursgal_lookup.csv.bz2"
    p = pytest._test_path / "data" / "BSA1_msgfplus_2021_03_22.mzid"
    db_path = pytest._test_path / "data" / "BSA.fasta"
    u = Unify(
        p,
        {
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
    assert isinstance(u.parser, MSGFPlus_2021_03_22_Parser)


def test_unify_get_msamanda_parser():
    rt_lookup_path = pytest._test_path / "data" / "BSA1_ursgal_lookup.csv.bz2"
    p = pytest._test_path / "data" / "BSA_msamanda_2_0_0_17442.csv"
    db_path = pytest._test_path / "data" / "BSA.fasta"
    u = Unify(
        p,
        {
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
    assert isinstance(u.parser, MSAmanda_2_Parser)


def test_unify_get_xtandem_parser():
    rt_lookup_path = pytest._test_path / "data" / "_ursgal_lookup.csv.bz2"
    p = pytest._test_path / "data" / "BSA1_xtandem_alanine.xml"
    db_path = pytest._test_path / "data" / "BSA.fasta"
    u = Unify(
        p,
        {
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
    assert isinstance(u.parser, XTandemAlanine_Parser)


def test_unify_get_msfragger_parser():
    rt_lookup_path = pytest._test_path / "data" / "_ursgal_lookup.csv.bz2"
    p = pytest._test_path / "data" / "test_Creinhardtii_QE_pH11_msfragger_3.tsv"
    db_path = pytest._test_path / "data" / "BSA.fasta"
    u = Unify(
        p,
        {
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
    assert isinstance(u.parser, MSFragger_3_Parser)


def test_unify_get_comet_parser():
    rt_lookup_path = pytest._test_path / "data" / "_ursgal_lookup.csv.bz2"
    p = pytest._test_path / "data" / "BSA1_comet_2020_01_4.mzid"
    db_path = pytest._test_path / "data" / "BSA.fasta"
    u = Unify(
        p,
        {
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
    assert isinstance(u.parser, Comet_2020_01_4_Parser)


def test_unify_get_mascot_parser():
    rt_lookup_path = pytest._test_path / "data" / "_ursgal_lookup.csv.bz2"
    p = pytest._test_path / "data" / "BSA1_mascot_2_6_2.dat"
    db_path = pytest._test_path / "data" / "BSA.fasta"
    u = Unify(
        p,
        {
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
    assert isinstance(u.parser, Mascot_2_6_2_Parser)
