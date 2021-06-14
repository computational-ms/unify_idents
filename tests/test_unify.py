#!/usr/bin/env python
from pathlib import Path
import uparma

from unify_idents.unify import Unify
from unify_idents.engine_parsers.omssa_parser import OmssaParser


def test_unify_determine_engine():
    rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
    p = Path(__file__).parent / "data" / "BSA1_mzml2mgf_0_0_1_omssa_2_1_9.csv_tmp"
    u = Unify(p, {"scan_rt_lookup_file": rt_lookup_path})
    assert u.engine == "omssa"


def test_unify_get_parser_classes():
    rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
    p = Path(__file__).parent / "data" / "BSA1_mzml2mgf_0_0_1_omssa_2_1_9.csv_tmp"
    u = Unify(p, {"scan_rt_lookup_file": rt_lookup_path})
    parsers = u._get_parser_classes()
    assert len(parsers) == 2  # currently we have dummy and omssa


def test_unify_get_omssa_parser():
    rt_lookup_path = Path(__file__).parent / "data" / "_ursgal_lookup.csv.bz2"
    p = Path(__file__).parent / "data" / "BSA1_mzml2mgf_0_0_1_omssa_2_1_9.csv_tmp"
    u = Unify(p, {"scan_rt_lookup_file": rt_lookup_path})
    parser = u._get_parser(p)
    assert isinstance(parser, OmssaParser)
