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
        self.fh = open(input_file)
        self.reader = iter(ElementTree.iterparse(self.fh, events=("end", "start")))
        self.peptide_lookup = self._get_peptide_lookup()

        self.cols_to_add = [
            "Raw data location",
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

    def __del__(self):
        self.fh.close()

    @classmethod
    def file_matches_parser(cls, file):
        ret_val = False
        p = Path(file)

        max_lines = 20

        if file.suffix == ".mzid":
            with open(file) as fin:
                mzml_iter = iter(ElementTree.iterparse(fin, events=("end", "start")))
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
        for n in self._next():
            u = self._unify_row(n)
            return u

    def _next(self):
        data = []
        while True:
            event, ele = next(self.reader, ("STOP", "STOP"))
            if event == "end" and ele.tag.endswith("SpectrumIdentificationResult"):
                for spec_result in list(
                    ele[::-1]
                ):  # iterate the list from end to start, since cvParams are after SpectrumIdentificationItem
                    if spec_result.tag.endswith("cvParam"):
                        if spec_result.attrib["name"] == "scan start time":
                            scan_time = spec_result.attrib["value"]
                        if spec_result.attrib["name"] == "scan number(s)":
                            spec_id = spec_result.attrib["value"]
                        if spec_result.attrib["name"] == "spectrum title":
                            spec_title = spec_result.attrib["value"]
                    else:
                        continue

                def all_items():
                    # todo use iterparse
                    for spec_result in ele:
                        if not spec_result.tag.endswith("SpectrumIdentificationItem"):
                            continue
                        pep_data = self.peptide_lookup[
                            spec_result.attrib["peptide_ref"]
                        ]
                        mods = []
                        for m in pep_data["Modifications"]:
                            name = m["name"]
                            pos = m["pos"]
                            mods.append(f"{name}:{pos}")
                        data = {
                            "Spectrum ID": spec_id,
                            "Retention Time (s)": scan_time,
                            "Peptide": pep_data["Sequence"],
                            "Modifications": ";".join(mods),
                            "Title": spec_title,
                            "Exp m/z": spec_result.attrib["experimentalMassToCharge"],
                            "Calc m/z": spec_result.attrib["calculatedMassToCharge"],
                            "Charge": spec_result.attrib["chargeState"],
                        }
                        for child in list(spec_result):
                            if child.tag.endswith("Param"):
                                n = child.attrib["name"]
                                if not n.startswith("MS-GF:"):
                                    n = f"MS-GF:{n}"
                                data[n] = child.attrib["value"]
                        yield data

                return all_items()
            if event == "STOP":
                raise StopIteration

    def _unify_row(self, row):
        for col_to_add in self.cols_to_add:
            if col_to_add not in row:
                row[col_to_add] = ""

        col_mapping = self.get_column_names(self.style)
        for new_key, old_key in col_mapping.items():
            if old_key in row.keys() and old_key != new_key:
                row[new_key] = row[old_key]
                del row[old_key]
        row["Search Engine"] = "MSGFPlus_2021_03_22"
        new_row = self.general_fixes(row)
        return UnifiedRow(**new_row)

    def _get_peptide_lookup(self):
        lookup = {}
        while True:
            event, ele = next(self.reader, ("STOP", "STOP"))
            if event == "end" and ele.tag.endswith("Peptide"):
                _id = ele.attrib.get("id", "")
                lookup[_id] = {}
                lookup[_id]["Modifications"] = []
                for child in ele:
                    if child.tag.endswith("PeptideSequence"):
                        seq = child.text
                        lookup[_id]["Sequence"] = seq
                    if child.tag.endswith("Modification"):
                        pos = child.attrib["location"]
                        mass = child.attrib["monoisotopicMassDelta"]
                        assert len(list(child)) == 1
                        name = list(child)[0].attrib["name"]
                        lookup[_id]["Modifications"].append(
                            {"pos": pos, "mass": mass, "name": name}
                        )
            if event == "start" and ele.tag.endswith("SpectrumIdentificationList"):
                break
        return lookup
