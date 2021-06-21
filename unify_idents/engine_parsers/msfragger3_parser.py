from unify_idents import UnifiedRow
from unify_idents.engine_parsers.base_parser import __BaseParser
from pathlib import Path
import csv


class MSFragger3Parser(__BaseParser):
    def __init__(self, input_file, params=None):
        super().__init__(input_file, params)
        if params is None:
            params = {}
        self.params = params
        self.input_file = input_file

        try:
            self.reader = csv.DictReader(open(input_file), delimiter="\t")
        except:
            self.reader = iter([])

        self.style = "msfragger_style_3"
        self.column_mapping = self.get_column_names()

        self.cols_to_add = [
            "Raw file location",
            "Spectrum Title",
            "uCalc m/z",
            "uCalc Mass",
            "Retention Time (s)",
            "Accuracy (ppm)",
            "Mass Difference",
            "Protein ID",
            "Sequence Start",
            "Sequence Stop",
            "Sequence Pre AA",
            "Sequence Post AA",
            "Enzyme Specificity",
            "Complies search criteria",
            "Conflicting uparam",
            "Search Engine",
        ]

    def __next__(self):
        line = next(self.reader)
        line = self._unify_row(line)
        return line

    def file_matches_parser(self):
        with open(self.input_file) as fin:
            headers = fin.readline()
            if set(headers.split()) == set(self.column_mapping.values()):
                ret_val = True
            else:
                ret_val = False
        return ret_val

    def _unify_row(self, row):
        new_row = {}
        for unify_name, omssa_name in self.column_mapping.items():
            new_row[unify_name] = row[omssa_name]
        for col in self.cols_to_remove:
            del new_row[col]
        for col in self.cols_to_add:
            if col not in new_row:
                print(f"ADD {col}")
                new_row[col] = ""
        new_row["Search Engine"] = "msfragger_3_0"
        new_row["Spectrum Title"] = "{file}.{specid}.{specid}.{charge}".format(
            file=self.params["Raw file location"],
            specid=new_row["Spectrum ID"],
            charge=new_row["Charge"],
        )
        new_row["Raw file location"] = self.params["Raw file location"]
        # breakpoint()
        new_row = self.general_fixes(new_row)
        return new_row

    def get_column_names(self):
        headers = self.param_mapper.get_default_params(style=self.style)[
            "header_translations"
        ]["translated_value"]
        return headers
