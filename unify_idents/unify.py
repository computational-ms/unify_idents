#!/usr/bin/env python
import bz2
import csv
import inspect
from peptide_mapper.mapper import UPeptideMapper
from chemical_composition import ChemicalComposition

from importlib import import_module
from pathlib import Path

import pandas as pd

PROTON = 1.00727646677


# inherit from pd.DataFrame
class UnifiedDataFrame:
    def __init__(self, rows, db_path, params=None):
        if params is None:
            self.params = {}
        else:
            self.params = params
        self.DELIMITER = "<|>"
        self.rows = rows
        self.df = pd.DataFrame(rows)
        # where do we get the database
        self.mapper = UPeptideMapper(db_path)
        self.cc = ChemicalComposition()
        all_peptides = [row["Sequence"] for row in rows]
        mapped_peptides = self.mapper.map_peptides(all_peptides)
        self.update_protein_mapping(mapped_peptides)

    def __len__(self):
        return len(self.df)

    def __iter__(self):
        return iter(self.rows)

    def update_protein_mapping(self, mapped_peptides):
        for _id, row in self.df.iterrows():
            mapping = mapped_peptides[row["Sequence"]]
            starts = []
            ids = []
            stops = []
            pre = []
            post = []
            for match in mapping:
                ids.append(match["id"])
                starts.append(str(match["start"]))
                stops.append(str(match["end"]))
                pre.append(str(match["pre"]))
                post.append(str(match["post"]))

            self.df.at[_id, "Is decoy"] = False
            for prot in ids:
                if self.params.get("decoy_tag", "decoy_") in prot:
                    self.df.at[_id, "Is decoy"] = True
                    break

            self.df.at[_id, "Protein ID"] = self.DELIMITER.join(ids)
            self.df.at[_id, "Sequence Pre AA"] = self.DELIMITER.join(pre)
            self.df.at[_id, "Sequence Post AA"] = self.DELIMITER.join(post)
            self.df.at[_id, "Sequence Start"] = self.DELIMITER.join(starts)
            self.df.at[_id, "Sequence Stop"] = self.DELIMITER.join(stops)

            # also do the mass calc here?
            self.cc.use(sequence=row["Sequence"], modifications=row["Modifications"])
            calc_mass = self.cc.mass()
            exp_mass = self.calc_mass(float(row["Exp m/z"]), int(row["Charge"]))
            charge = int(row["Charge"])
            calc_mz = self.calc_mz(calc_mass, charge)
            self.df.at[_id, "uCalc m/z"] = calc_mz
            self.df.at[_id, "uCalc Mass"] = calc_mass + PROTON
            if self.df.at[_id, "Calc m/z"] == "":
                self.df.at[_id, "Calc m/z"] = calc_mz
            acc = (
                (float(row["Exp m/z"]) - float(self.df.at[_id, "uCalc m/z"]))
                / self.df.at[_id, "uCalc m/z"]
                * 1e6
            )
            self.df.at[_id, "Accuracy (ppm)"] = acc

        return

    def calc_mz(self, mass, charge):
        return (mass + (charge * PROTON)) / charge

    def calc_mass(self, mz, charge):
        calc_mass = mz * charge - (charge * PROTON)
        return calc_mass


# inherit from ordered dict?
class UnifiedRow:
    def __init__(self, **kwargs):
        self.columns = kwargs.keys()
        self.data = kwargs
        self.col_order = [
            "Spectrum Title",
            "Spectrum ID",
            "Sequence",
            "Modifications",
            "Charge",
            "Protein ID",
            "Retention Time (s)",
            "Exp m/z",
            "Calc m/z",
            "uCalc m/z",
            "uCalc Mass",
            "Accuracy (ppm)",
            "Mass Difference",
            "Sequence Start",
            "Sequence Stop",
            "Sequence Pre AA",
            "Sequence Post AA",
            "Enzyme Specificity",
            "Complies search criteria",
            "Conflicting uparma",
            "Search Engine",
        ]
        self.string_repr = None

    def __getitem__(self, key):
        return self.data[key]

    def __setitem__(self, key, value):
        if key not in self.data:
            raise KeyError("Cant add new key")
        self.data[key] = value

    def __contains__(self, key):
        return key in self.data

    def __str__(self):
        # needs fix, only return unify cols and not Engine:name columns
        if self.string_repr is None:
            self.string_repr = []
            for col in self.col_order:
                self.string_repr.append(str(self.data.get(col, "")))
            self.string_repr = ", ".join(self.string_repr)
        return self.string_repr

    def to_dict(self):
        return self.data


class Unify:
    def __init__(self, input_file, params=None):
        """Summary"""
        self.input_file = input_file
        self.params = params
        if not isinstance(self.input_file, Path):
            self.input_file = Path(self.input_file)
        else:
            self.input_file = self.input_file
        if params is None:
            self.params = {}
        else:
            self.params = params

        self.parser = self._get_parser(self.input_file)

    def __iter__(self):
        return self

    def __next__(self):
        line = next(self.parser)
        return line

    def _get_parser_classes(self):
        classes = []
        for eng_type in ["quant", "ident"]:
            fstring = "unify_idents.engine_parsers.{eng_type}.{module}"
            all_modules = []
            p = Path(__file__).parent / "engine_parsers" / eng_type
            for child in p.iterdir():
                if str(child).endswith(".py") and not str(child.stem).startswith("__"):
                    name = child.stem
                    formatted_string = fstring.format(module=name, eng_type=eng_type)
                    members = inspect.getmembers(import_module(formatted_string))
                    for name, obj in dict(members).items():
                        if (
                            inspect.isclass(obj)
                            and obj.__module__ == formatted_string
                            and not "__" in name
                        ):
                            classes.append(obj)
        return classes

    def _get_parser(self, input_file):
        all_parsers = self._get_parser_classes()
        for parser_class in all_parsers:
            if parser_class.file_matches_parser(input_file) is True:
                break
        parser = parser_class(input_file, params=self.params)
        return parser

    def read_rt_lookup_file(self, scan_rt_lookup_path):
        with bz2.open(scan_rt_lookup_path, "rt") as fin:
            lookup = {}
            reader = csv.DictReader(fin)
            for line in reader:
                lookup.setdefault(
                    line["File"], {"scan2rt": {}, "rt2scan": {}, "scan2mz": {}}
                )
                file, scan, rt, mz = (
                    line["File"],
                    line["Spectrum ID"],
                    line["RT"],
                    line["Precursor mz"],
                )
                lookup[file]["scan2rt"][int(scan)] = float(rt)
                lookup[file]["rt2scan"][float(rt)] = int(scan)
                lookup[file]["scan2mz"][int(scan)] = float(mz)
        return lookup

    def get_dataframe(self):
        data = []
        for unified_row in self:
            if unified_row is not None:
                data.append(unified_row.to_dict())
        return UnifiedDataFrame(data, db_path=self.params["database"])
