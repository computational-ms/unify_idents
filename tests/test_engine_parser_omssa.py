#!/usr/bin/env python
from pathlib import Path
from unify_idents.engine_parsers.omssa_parser import OmssaParser
import uparma


def test_engine_parsers_omssa_init():
    input_file = (
        Path(__file__).parent / "data" / "BSA1_mzml2mgf_0_0_1_omssa_2_1_9.csv_tmp"
    )
    parser = OmssaParser(input_file)


def test_engine_parsers_omssa_unify_row():
    input_file = (
        Path(__file__).parent / "data" / "BSA1_mzml2mgf_0_0_1_omssa_2_1_9.csv_tmp"
    )
    parser = OmssaParser(input_file)
    for row in parser:
        print(row)
    assert 1 == 2
