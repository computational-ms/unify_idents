#!/usr/bin/env python
import csv

import uparma

from unify_idents import UnifiedRow
from unify_idents.engine_parsers.base_parser import __BaseParser

import xml.etree.ElementTree as ElementTree
from xml.etree.ElementTree import ParseError
from pathlib import Path


class MSGFPlus_2021_03_22(__BaseParser):
    def __init__(self, input_file, params=None):
        super().__init__(input_file, params)
        if params is None:
            params = {}
        self.params = params
        self.input_file = input_file

        self.style = "msgfplus_style_1"
        self.column_mapping = self.get_column_names(self.style)
        fh = open(input_file)
        self.reader = iter(ElementTree.iterparse(fh, events=("end",)))

    def __del__(self):
        fh.close()

    @classmethod
    def file_matches_parser(cls, file):
        ret_val = False
        p = Path(file)

        max_lines = 20
        with open(file) as fin:
            mzml_iter = iter(ElementTree.iterparse(fin, events=("end",)))
            for pos, (event, ele) in enumerate(mzml_iter):
                if pos > max_lines:
                    ret_val = False
                    break
                if ele.tag.endswith("AnalysisSoftware"):
                    name = ele.attrib.get("name", "")
                    version = ele.attrib.get("version", "")
                    if name == "MS-GF+" and version == "Release (v2021.03.22)":
                        ret_val = True
                        break
        return ret_val

    def get_column_names(self, style):
        headers = self.param_mapper.get_default_params(style=style)[
            "header_translations"
        ]["translated_value"]
        return headers

    def __iter__(self):
        return self

    def __next__(self):
        n = next(self.reader)
        u = self._unify_row(n)
        return u

    def _unify_row(self, row):
        pass
