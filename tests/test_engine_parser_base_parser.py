from pathlib import Path

import pytest

from unify_idents.engine_parsers.base_parser import BaseParser


def test_base_parser_read_rt_lookup_file():
    rt_lookup_path = Path(__file__).parent / "data" / "BSA1_ursgal_lookup.csv.bz2"
    input_file = (
        Path(__file__).parent / "data" / "test_Creinhardtii_QE_pH11_xtandem_alanine.xml"
    )

    bp = BaseParser(input_file, params={"rt_pickle_name": rt_lookup_path})
    rt_lookup = bp._read_rt_lookup_file()
    assert (rt_lookup["Unit"] == 1).all()
    assert pytest.approx(
        rt_lookup["Precursor mz"].mean(), 550.8444810049874
    )  # check consistency
    assert pytest.approx(rt_lookup.loc[2450, "Precursor mz"], 618.2697)
    assert pytest.approx(rt_lookup.loc[2450, "RT"], 92067.714)
