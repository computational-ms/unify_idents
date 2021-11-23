"""Quant parser."""
import pandas as pd
import re
from pathlib import Path
from unify_idents.engine_parsers.base_parser import QuantBaseParser


class FlashLFQ_1_2_0_Parser(QuantBaseParser):
    """File parser for Flash LFQ."""

    def __init__(self, *args, **kwargs):
        """Initialize parser.

        Reads in data file and provides mappings.
        """
        super().__init__(*args, **kwargs)
        self.style = "flash_lfq_style_1"
        self.mapping_dict = {
            v: k
            for k, v in self.param_mapper.get_default_params(style=self.style)[
                "header_translations"
            ]["translated_value"].items()
        }
        self.df = pd.read_csv(self.input_file, delimiter="\t")
        self.df.rename(columns=self.mapping_dict, inplace=True)

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
        # raise NotImplementedError
        # do column conversion here
        # breakpoint()
        self.df["Spectrum ID"] = -1
        self.df["Linked Spectrum ID"] = -1
        self.df["Raw Quant Value"] = -1
        self.df["FWHM"] = ""
        self.df["Label"] = "LabelFree"
        self.df["Condition"] = self.df["Raw Filename"]
        self.df["Raw Filename"] = self.df["Raw Filename"].map(
            lambda x: Path(self.params.get("Raw data location", "")) / Path(x).stem
        )
        self.df["Quant Group"] = ""
        self.df["Score"] = ""
        self.df["Processing Level"] = "ChromatographicPeak"
        self.df["Quant Run ID"] = "FlashLFQ"
        self.df["Coalescence"] = ""
        self.process_unify_style()
        return self.df

    def translate_mods(self, full_sequence):
        """Extract modifications from full_sequence and format as {mod_1}:{pos_1};{mod_n}:{pos_n}.

        Args:
            full_sequence (str): sequence including mod (e.g. ELVISC[Carbamidomethyl]M[Oxidation])

        Returns:
            str: extracted mods as described in summary
        """
        # TODO extract C-terminal mods, position should be seq_len + 1
        cumulative_match_length = 0
        regex = re.compile("\[(.*?)\]")
        mods = []
        for match in regex.finditer(full_sequence):
            mods.append(f"{match.group(1)}:{match.start() - cumulative_match_length}")
            cumulative_match_length += len(match.group())
        return ";".join(mods)
