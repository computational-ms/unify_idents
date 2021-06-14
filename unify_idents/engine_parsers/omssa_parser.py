#!/usr/bin/env python
import csv
from unify_idents import (
    UnifiedRow,
)
import uparma
from unify_idents.engine_parsers.base_parser import __BaseParser


"""OMSSA
        * Carbamidomethyl is updated and set
"""
# add a base parser, which has all the common functionality?
# base parser could also do type casting
# use ABC?
# sanity check
class OmssaParser(__BaseParser):
    def __init__(self, input_file, params=None):
        super().__init__(input_file, params)
        if params is None:
            params = {}
        self.params = params
        self.input_file = input_file

        self.reader = csv.DictReader(open(input_file))
        self.style = "omssa_style_1"
        self.column_mapping = self.get_column_names()

        # copy pasta from ursgal1 unify_csv
        self.cols_to_remove = [
            "proteinacc_start_stop_pre_post_;",
            "Start",
            "Stop",
            "NIST score",
            "gi",
            "Accession",
        ]

        self.cols_to_add = [
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

    def file_matches_parser(self):
        # TODO implement file sensing
        return True

    def __iter__(self):
        return self

    def __next__(self):
        n = next(self.reader)
        u = self._unify_row(n)
        return u

    def _unify_row(self, row):
        new_row = {}
        for unify_name, omssa_name in self.column_mapping.items():
            new_row[unify_name] = row[omssa_name]
        for col in self.cols_to_remove:
            del new_row[col]
        for col in self.cols_to_add:
            new_row[col] = ""
        # breakpoint()
        new_row["Spectrum ID"] = int(new_row["Spectrum Title"].split(".")[1])
        new_row["Search Engine"] = "omssa_2_1_9"
        new_row = self.general_fixes(new_row)
        breakpoint()
        # TODO add fixed modifications
        # TODO
        # all the other transformations
        return UnifiedRow(**new_row)

    def get_column_names(self):
        # create own uparma mapper
        headers = self.param_mapper.get_default_params(style=self.style)[
            "header_translations"
        ]["translated_value"]
        return headers
