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


col_mapping = {"seq": "Sequence", "z": "Charge", "hyperscore": "X!Tandem:Hyperscore"}


col_mapping = {
    "X!Tandem:delta": "delta",
    "X!Tandem:nextscore": "nextscore",
    "X!Tandem:y_score": "y_score",
    "X!Tandem:y_ions": "y_ions",
    "X!Tandem:b_score": "b_score",
    "X!Tandem:b_ions": "b_ions",
    "Sequence": "seq",
    "Charge": "z",
    "X!Tandem:Hyperscore": "hyperscore",
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

        self.cols_to_remove = [
            "id",
            "start",
            "end",
            "pre",
            "post",
            "missed_cleavages",
            "expect",
        ]

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
                        row["Spectrum Title"] = spec_title.split()[0]
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
        for new_col, old_col in col_mapping.items():
            row[new_col] = row[old_col]
            del row[old_col]
        for col in self.cols_to_remove:
            del row[col]
        modstring = self.map_mod_names(row)

        row["Modifications"] = modstring
        row["Search Engine"] = "xtandem_alanine"
        row["Spectrum ID"] = row["Spectrum Title"].split(".")[1]
        row = self.general_fixes(row)
        return UnifiedRow(**row)
