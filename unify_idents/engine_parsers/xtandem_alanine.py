from unify_idents import UnifiedRow
from unify_idents.engine_parsers.base_parser import __BaseParser
from pathlib import Path
import re
import csv
from decimal import Decimal, getcontext, ROUND_UP
import itertools
from loguru import logger
import xml.etree.ElementTree as ElementTree
import copy


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
        while True:
            try:
                gen = self._next()
                for x in gen():
                    x = self._unify_row(x)
                    yield x
            except StopIteration:
                break

    def __next__(self):
        return next(self.__iter__())

    def _next(self):
        while True:
            event, element = next(self.xml_iter, ("STOP", "STOP"))
            if event == "STOP":
                raise StopIteration
            if (
                event == "start"
                and element.tag.endswith("group")
                and "z" in element.attrib
            ):
                charge = element.attrib["z"]
                prec_mz = element.attrib["mh"]

            if (
                event == "end"
                and element.tag.endswith("group")
                and "z" in element.attrib
                # and element.attrib["id"] == "552"
            ):
                spec_title = element.findall('.//**[@label="Description"]')[0].text
                charge = element.attrib["z"]

                def result_iterator():
                    for child in element.findall(".//protein"):
                        # which mh is Exp m/z and which is Calc m/z??
                        domain = child.findall(".//domain")[0]
                        row = copy.copy(domain.attrib)
                        row["Exp m/z"] = prec_mz
                        row["Calc m/z"] = row["mh"]
                        del row["mh"]
                        row["Modifications"] = []
                        row["Spectrum Title"] = spec_title
                        row["z"] = charge
                        mods = domain.findall(".//aa")
                        for m in mods:
                            mass, abs_pos = m.attrib["modified"], m.attrib["at"]
                            # abs pos is pos in protein, rel pos is pos in peptide
                            rel_pos = int(abs_pos) - int(row["start"])
                            row["Modifications"].append(f"{mass}:{rel_pos}")
                        yield row

                return result_iterator

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
