#!/usr/bin/env python
from pathlib import Path
import uparma

from unify_idents.unify import Unify
from unify_idents.engine_parsers.omssa_parser import OmssaParser


def test_unify_determine_engine():
    d = {
        "ufiles": (
            Path(__file__).parent / "data" / "BSA1_mzml2mgf_0_0_1_omssa_2_1_9.csv_tmp",
        )
    }
    u = Unify(**d)
    assert u.engine == "omssa"


def test_unify_get_parser_classes():
    d = {
        "ufiles": (
            Path(__file__).parent / "data" / "BSA1_mzml2mgf_0_0_1_omssa_2_1_9.csv_tmp",
        )
    }
    u = Unify(**d)
    parsers = u._get_parser_classes()
    assert len(parsers) == 2  # currently we have dummy and omssa


def test_unify_get_omssa_parser():
    d = {
        "ufiles": (
            Path(__file__).parent / "data" / "BSA1_mzml2mgf_0_0_1_omssa_2_1_9.csv_tmp",
        )
    }
    u = Unify(**d)
    parser = u._get_parser(d["ufiles"][0])
    assert isinstance(parser, OmssaParser)
