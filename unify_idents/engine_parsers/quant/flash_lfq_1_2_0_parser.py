import csv
import difflib
import re
from pathlib import Path

from loguru import logger

from unify_idents.engine_parsers.base_parser import __QuantBaseParser

col_mapping = {
    "File Name": "file_name",
    "Base Sequence": "trivial_name",
    # "Protein Group": "Protein IDs",
    # "Peptide Monoisotopic Mass": "Mass",
    "Peak RT Apex": "retention_time",
    "Precursor Charge": "charge",
    # "Theoretical MZ": "Calc m/z",
    "Peak intensity": "quant_value",
    "Peak Apex Mass Error (ppm)": "delta_mz",
}


class FlashLFQ(__QuantBaseParser):
    def __init__(self, input_file, params=None):
        """Summary

        Args:
            input_file (str): FlashLFQ QuantifiedPeaks.tsv
            params (None, optional): FlashLFQ specific parameters
        """
        super().__init__(input_file, params)
        self.style = "flash_lfq_style_1"
        self.column_mapping = self.get_column_names()
        try:
            self.reader = csv.DictReader(open(input_file), delimiter="\t")
        except IOError:
            self.reader = iter([])

    @classmethod
    def file_matches_parser(cls, file):
        """Check if `file` is a valid input file for this class

        Args:
            file (str): input file path

        Returns:
            bool: Wether the input file is a valida input or not
        """
        flash_lfq_columns = set(
            [
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
            ]
        )
        ret_val = False
        with open(file) as fin:
            file_column_names = set(fin.readline().strip().split("\t"))
        if flash_lfq_columns == file_column_names:
            ret_val = True
        return ret_val

    def __iter__(self):
        return self

    def __next__(self):
        line = next(self.reader)
        line = self._unify_row(line)
        return line

    def _unify_row(self, row):
        """Transform row into unified_quant format.

        Args:
            row (dict): Row as present in input file

        Returns:
            dict: Transformed row
        """
        self.cc.use(
            sequence=row["Base Sequence"],
            modifications=self.extract_mods(row["Full Sequence"]),
        )  # currently only works for peptides
        new_row = {}
        for flash_name, unify_name in col_mapping.items():
            new_row[unify_name] = row[flash_name]

        keys = list(row.keys())
        for old_key in keys:
            if old_key not in col_mapping.keys() and old_key != "":
                new_row[f"FlashLFQ:{old_key}"] = row[old_key]
            del row[old_key]
        calc_mz = self.calc_mz(self.cc.mass(), int(new_row["charge"]))

        # Identifier
        new_row["spectrum_id"] = ""
        new_row["chemical_composition"] = self.cc.hill_notation_unimod()

        # spectrum meta data
        new_row["precursor_spectrum_id"] = ""
        new_row["ms_level"] = 1

        # Quant data
        new_row["quant_engine"] = "FlashLFQ_1_2_0"
        new_row["quant_score"] = ""
        new_row["quant_group"] = ""
        new_row["processing_level"] = "gaussian_apex"

        # Quant meta data
        new_row["label"] = "LabelFree"
        new_row["condition"] = ""
        new_row["ident_reference"] = ""

        # Derived data
        new_row["fwhm"] = ""
        new_row["s2i"] = ""
        new_row["p2t"] = ""
        new_row["coalescence"] = ""

        if (raw_name := self.params.get("Raw file location", None)) is not None:
            new_row["file_name"] = raw_name

        if self.check_required_headers(new_row) is False:
            logger.error("Not all required headers are present!")
        return new_row

    def extract_mods(self, full_sequence):
        """Extract modifications from full_sequence and format as {mod_1}:{pos_1};{mod_n}:{pos_n}

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
