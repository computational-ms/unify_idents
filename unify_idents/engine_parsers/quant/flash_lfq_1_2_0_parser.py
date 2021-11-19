"""Quant parser."""
import pandas as pd

from unify_idents.engine_parsers.base_parser import __QuantBaseParser


class FlashLFQ_1_2_0_Parser(__QuantBaseParser):
    """File parser for Flash LFQ."""

    def __init__(self, *args, **kwargs):
        """Initialize parser.

        Reads in data file and provides mappings.
        """
        super().__init__(*args, **kwargs)
        self.style = "flash_lfq_style_1"
        self.reference_columns = self.param_mapper.get_default_params(style=self.style)[
            "header_translations"
        ]["translated_value"]
        self.df = pd.read_csv(self.input_file, delimiter="\t")

    @classmethod
    def check_parser_compatibility(cls, file):
        """Assert compatibility between file and parser.

        Args:
            file (str): path to input file

        Returns:
            bool: True if parser and file are compatible

        """
        is_tsv = file.as_posix().endswith(".tsv")
        flash_lfq_columns = {
            "File Name",
            "Base Sequence",
            "Full Sequence",
            "Protein Group",
            "Peptide Monoisotopic Mass",
            "MS2 Retention Time",
            "Precursor Charge",
            "Theoretical MZ",
            "Peak intensity",
            "Peak RT Start",
            "Peak RT Apex",
            "Peak RT End",
            "Peak MZ",
            "Peak Charge",
            "Num Charge States Observed",
            "Peak Detection Type",
            "MBR Score",
            "PSMs Mapped",
            "Base Sequences Mapped",
            "Full Sequences Mapped",
            "Peak Split Valley RT",
            "Peak Apex Mass Error (ppm)",
        }
        with open(file.as_posix()) as f:
            head = set(f.readline().replace("\n", "").split("\t"))
        headers_match = len(flash_lfq_columns.difference(head)) == 0
        return is_tsv and headers_match

    def unify(self):
        """Primary method to read and unify engine output.

        Returns:
            self.df (pd.DataFrame): unified dataframe
        """
        # TODO: fill this
        raise NotImplementedError

        return self.df
