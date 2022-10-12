"""Quant base parser class."""

from chemical_composition import ChemicalComposition

from unify_idents.engine_parsers.base_parser import BaseParser


class QuantBaseParser(BaseParser):
    """Base class of all quant parsers."""

    def __init__(self, *args, **kwargs):
        """Initialize parser.

        Reads in data file and provides mappings.
        """
        super().__init__(*args, **kwargs)
        self.cc = ChemicalComposition()
        self.required_headers = {
            "file_name",
            "spectrum_id",
            "trivial_name",
            "chemical_composition",
            "precursor_spectrum_id",
            "retention_time",
            "charge",
            "quant_run_id",
            "quant_value",
            "quant_score",
            "quant_group",
            "processing_level",
            "delta_mz",
            "label",
            "condition",
            "ident_reference",
            "fwhm",
            "s2i",
            "p2t",
            "coalescence",
        }

    def process_unify_style(self):
        """Apply sanitizing methods."""
        pass