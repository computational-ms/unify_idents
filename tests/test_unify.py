#!/usr/bin/env python
from pathlib import Path

import pandas as pd
import pytest

from unify_idents.engine_parsers.ident.msamanda_2_parser import MSAmanda_2_Parser
from unify_idents.engine_parsers.ident.msgfplus_2021_03_22_parser import (
    MSGFPlus_2021_03_22,
)
from unify_idents.engine_parsers.ident.omssa_2_1_9_parser import OmssaParser
from unify_idents.engine_parsers.ident.xtandem_alanine import XTandemAlanine
from unify_idents.engine_parsers.ident.msfragger_3_parser import MSFragger3Parser
from unify_idents.engine_parsers.ident.comet_2020_01_4_parser import (
    Comet_2020_01_4_Parser,
)
from unify_idents.unify import Unify


def test_unify_get_parser_classes():
    rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
    p = Path(__file__).parent / "data" / "test_Creinhardtii_QE_pH11_omssa_2_1_9.csv"
    db_path = Path(__file__).parent / "data" / "BSA.fasta"
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
            "omssa_mod_dir": Path(__file__).parent / "data",
        },
    )
    # currently msamanda, msfragger, msgfplus, omssa, xtandem, comet, flash_lfq, and dummy
    assert len(u._parser_classes) == 8


def test_unify_get_omssa_parser():
    rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
    p = Path(__file__).parent / "data" / "test_Creinhardtii_QE_pH11_omssa_2_1_9.csv"
    db_path = Path(__file__).parent / "data" / "BSA.fasta"
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
            "omssa_mod_dir": Path(__file__).parent / "data",
        },
    )
    assert isinstance(u.parser, OmssaParser)


def test_unify_get_msgfplus_parser():
    rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
    p = Path(__file__).parent / "data" / "BSA1_msgfplus_2021_03_22.mzid"
    db_path = Path(__file__).parent / "data" / "BSA.fasta"
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
    assert isinstance(u.parser, MSGFPlus_2021_03_22)


def test_unify_get_msamanda_parser():
    rt_lookup_path = Path(__file__).parent / "data" / "BSA1_ursgal_lookup.csv.bz2"
    p = Path(__file__).parent / "data" / "BSA_msamanda_2_0_0_17442.csv"
    db_path = Path(__file__).parent / "data" / "BSA.fasta"
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
    rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
    p = Path(__file__).parent / "data" / "BSA1_xtandem_alanine.xml"
    db_path = Path(__file__).parent / "data" / "BSA.fasta"
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
    assert isinstance(u.parser, XTandemAlanine)


def test_unify_get_msfragger_parser():
    rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
    p = Path(__file__).parent / "data" / "test_Creinhardtii_QE_pH11_msfragger_3.tsv"
    db_path = Path(__file__).parent / "data" / "BSA.fasta"
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
    assert isinstance(u.parser, MSFragger3Parser)


def test_unify_get_comet_parser():
    rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
    p = Path(__file__).parent / "data" / "BSA1_comet_2020_01_4.mzid"
    db_path = Path(__file__).parent / "data" / "BSA.fasta"
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
    pass
