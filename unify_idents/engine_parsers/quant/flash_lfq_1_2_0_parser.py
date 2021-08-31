import csv
import difflib
import re
from pathlib import Path

from loguru import logger

from unify_idents.engine_parsers.base_parser import __QuantBaseParser

col_mapping = {
    "File Name": "Raw data location",
    "Base Sequence": "Sequence",
    "Protein Group": "Protein IDs",
    "Peptide Monoisotopic Mass": "Mass",
    "Peak RT Apex": "Retention Time (s)",
    "Precursor Charge": "Charge",
    "Theoretical MZ": "Calc m/z",
    "Peak intensity": "Quant Value",
    "Peak Apex Mass Error (ppm)": "PPM",
}


class FlashLFQ(__QuantBaseParser):
    def __init__(self, input_file, params=None):
        super().__init__(input_file, params)
        self.style = "flash_lfq_style_1"
        self.column_mapping = self.get_column_names()
        try:
            self.reader = csv.DictReader(open(input_file), delimiter="\t")
        except IOError:
            self.reader = iter([])

    @classmethod
    def file_matches_parser(cls, file):
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
        new_row = {}
        for flash_name, unify_name in col_mapping.items():
            new_row[unify_name] = row[flash_name]
        self.cc.use(
            sequence=new_row["Sequence"],
            modifications=self.extract_mods(row["Full Sequence"]),
        )  # currently only works for peptides
        calc_mz = self.calc_mz(float(new_row["Mass"]), int(new_row["Charge"]))
        new_row["Spectrum ID"] = -1
        new_row["Linked Spectrum ID"] = -1
        new_row["Chemical Composition"] = self.cc.hill_notation_unimod()
        new_row["Raw Quant Value"] = -1
        new_row["MZ Delta"] = float(new_row["PPM"]) * 1e-6 * calc_mz
        new_row["FWHM"] = ""
        new_row["Label"] = "LabelFree"
        new_row["Condition"] = Path(new_row["Raw data location"]).stem
        new_row["Quant Group"] = ""
        new_row["Score"] = ""
        new_row["Processing Level"] = "ChromatographicPeak"
        new_row["Quant Run ID"] = "XXX+FlashLFQ"
        new_row["Coalescence"] = ""
        # breakpoint()
        if self.check_required_headers(new_row) is False:
            logger.error("Not all required headers are present!")
        return new_row

    def extract_mods(self, full_sequence):
        cumulative_match_length = 0
        regex = re.compile("\[(.*?)\]")
        mods = []
        for match in regex.finditer(full_sequence):
            mods.append(f"{match.group(1)}:{match.start() - cumulative_match_length}")
            cumulative_match_length += len(match.group())
        return ";".join(mods)
