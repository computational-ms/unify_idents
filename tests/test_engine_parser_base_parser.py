import pytest

from unify_idents.engine_parsers.base_parser import BaseParser


def test_uninitialized_parser_compatiblity_is_false():
    input_file = (
        pytest._test_path / "data" / "test_Creinhardtii_QE_pH11_xtandem_alanine.xml"
    )
    compat = BaseParser.check_parser_compatibility(input_file)
    assert compat is False
