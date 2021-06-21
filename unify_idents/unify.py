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
    def __init__(self, rows, db_path):
        self.DELIMITER = "<|>"
        self.rows = rows
        self.df = pd.DataFrame(rows)
        # where do we get the database
        self.mapper = UPeptideMapper(db_path)
        self.cc = ChemicalComposition()
        all_peptides = [row["Sequence"] for row in rows]
        mapped_peptides = self.mapper.map_peptides(all_peptides)
        self.update_protein_mapping(mapped_peptides)

    def __iter__(self):
        return iter(self.rows)

    def update_protein_mapping(self, mapped_peptides):
        for _id, row in self.df.iterrows():
            mapping = mapped_peptides[row["Sequence"]]
            # breakpoint()
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

            self.df.at[_id, "Protein ID"] = self.DELIMITER.join(ids)
            self.df.at[_id, "Sequence Pre AA"] = self.DELIMITER.join(pre)
            self.df.at[_id, "Sequence Post AA"] = self.DELIMITER.join(post)
            self.df.at[_id, "Sequence Start"] = self.DELIMITER.join(starts)
            self.df.at[_id, "Sequence Stop"] = self.DELIMITER.join(stops)

            # also do the mass calc here?
            self.cc.use(sequence=row["Sequence"], modifications=row["Modifications"])
            mass = self.cc.mass()
            charge = int(row["Charge"])
            mz = self.calc_mz(mass, charge)
            # exp m/z is actually the math
            engine_mass = float(row["Exp m/z"])
            self.df.at[_id, "uCalc m/z"] = mz
            self.df.at[_id, "uCalc Mass"] = mass
            self.df.at[_id, "Mass Difference"] = engine_mass - mass
            self.df.at[_id, "Accuracy (ppm)"] = (engine_mass - mass) / mass / 5e-6
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

    # def __repr__(self):
    #     return self.__str__()

    # def calc_mz(self):
    #     # use chemical composition
    #     # do we really want a new CC object for every row?
    #     self.data["uCalc m/z"] = 0
    #     self.data["uCalc mass"] = 0
    #     self.data["Mass Difference"] = 0
    #     self.data["Accuracy (ppm)"] = 0


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

        self.scan_rt_path = self.params.get("scan_rt_lookup_file", None)
        if self.scan_rt_path is None:
            raise Exception("Meaningfull error message")
        else:
            self.scan_rt_lookup = self.read_rt_lookup_file(self.scan_rt_path)
        self.parser = self._get_parser(self.input_file)

    def __iter__(self):
        return self

    def __next__(self):
        line = next(self.parser)
        # TODO
        # line = self.parser.general_fixes(line)
        # do some magic here, like calling methods of row (e.g. calc_mz)
        return line

    def _get_parser_classes(self):
        classes = []
        fstring = "unify_idents.engine_parsers.{module}"
        all_modules = []
        p = Path(__file__).parent / "engine_parsers"
        for child in p.iterdir():
            if str(child).endswith(".py") and not str(child.stem).startswith("__"):
                name = child.stem
                formatted_string = fstring.format(module=name)
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
            parser = parser_class(input_file, params=self.params)
            if parser.file_matches_parser() is True:
                break
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
            data.append(unified_row.to_dict())
        return UnifiedDataFrame(data, db_path=self.params["database"])
