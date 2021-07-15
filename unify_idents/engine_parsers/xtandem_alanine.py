from unify_idents import UnifiedRow
from unify_idents.engine_parsers.base_parser import __BaseParser
from pathlib import Path
import re
import csv
from decimal import Decimal, getcontext, ROUND_UP
import itertools
from loguru import logger
import xml.etree.ElementTree as ElementTree


col_mapping = {
    "seq": "Sequence",
    "z": "Charge",
}


class XTandemAlanine(__BaseParser):
    def __init__(self, input_file, params=None):
        super().__init__(input_file, params)
        if params is None:
            params = {}
        self.params = params
        self.input_file = input_file
        self.fin = open(input_file)
        self.xml_iter = iter(ElementTree.iterparse(self.fin, events=("end", "start")))
        self.raw_file = None

    def __del__(self):
        try:
            self.fin.close()
        except NameError:
            pass

    @classmethod
    def file_matches_parser(cls, file):
        ret_val = False

        if not str(file).endswith(".xml"):
            ret_val = False
        else:
            with open(file) as fout:
                for i, line in enumerate(fout):
                    if i > 10:  # only check first ten lines
                        break
                    if (
                        # '<?xml-stylesheet type="text/xsl" href="tandem-style.xsl"?>'
                        "tandem-style.xsl"
                        in line
                    ):
                        ret_val = True
                        break
        return ret_val

    def __iter__(self):
        return self

    def __next__(self):
        line = self._next()
        unified_line = self._unify_row(line)
        return unified_line

    def _next(self):
        row = {"Modifications": set()}

        while True:
            event, element = next(self.xml_iter, ("STOP", "STOP"))
            if event == "STOP":
                raise StopIteration
            elif (
                event == "start"
                and element.tag.endswith("group")
                and element.attrib["type"] == "model"
            ):
                # breakpoint()
                element.attrib["Exp m/z"] = (
                    float(element.attrib["mh"]) / int(element.attrib["z"]) + self.PROTON
                )
                row.update(element.attrib)
            elif event == "start" and element.tag.endswith("domain"):
                element.attrib["Calc m/z"] = self.calc_mz(
                    float(element.attrib["mh"]), float(row["z"])
                )
                row.update(element.attrib)
            elif event == "end" and element.tag.endswith("aa"):
                mass = element.attrib["modified"]
                pos = int(element.attrib["at"]) - int(row["start"])
                row["Modifications"].add(f"{mass}:{pos}")
                # if row["seq"] == "MLNMLIVFRFLRIIPSMK":
                #     breakpoint()
            elif (
                element.tag.endswith("note")
                and element.attrib["label"] == "Description"
            ):
                row["Spectrum Title"] = element.text
            elif (
                event == "start"
                and element.tag.endswith("bioml")
                and element.attrib["label"].startswith("models from")
            ):
                self.raw_file = element.attrib["label"].split("'")[1]
                # breakpoint()
            elif (
                event == "end"
                and element.tag.endswith("group")
                and element.attrib["type"] == "model"
            ):
                row["Raw data location"] = self.raw_file
                row["Modifications"] = list(row["Modifications"])
                # breakpoint()
                return row

    def _unify_row(self, row):
        for old_col, new_col in col_mapping.items():
            row[new_col] = row[old_col]
            del row[old_col]
        modstring = self.map_mod_names(row)

        row["Modifications"] = modstring
        row["Search Engine"] = "X!TandemAlanine"
        row["Spectrum ID"] = row["Spectrum Title"].split(".")[1]
        row = self.general_fixes(row)
        return UnifiedRow(**row)
