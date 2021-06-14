#!/usr/bin/env python
import inspect
from importlib import import_module
from pathlib import Path
import bz2
import csv

# Add UnifiedTable


class UnifiedDataFrame:
    pass


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
            "uCalc_m/z",
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
        if self.string_repr is None:
            self.string_repr = []
            for col in self.col_order:
                self.string_repr.append(self.data.get(col, ""))
            self.string_repr = ", ".join(self.string_repr)
        return self.string_repr

    def __repr__(self):
        return self.__str__()

    def calc_mz(self):
        # use chemical composition
        pass


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
        # self.engine = self._determine_engine(self.input_file)
        self.parser = self._get_parser(self.input_file)

    def __iter__(self):
        return self

    def __next__(self):
        line = next(self.parser)
        # do some magic here, like calling methods of row (e.g. calc_mz)
        return line

    # def _determine_engine(self, input_file):
    #     """Determine the engine used for generating `input_file`

    #     Args:
    #         input_file (str): path to the input file
    #     """
    #     engine = None
    #     if input_file.suffix == ".xml":
    #         # either xtandem or msgfplus
    #         pass
    #     elif input_file.suffix == ".csv_tmp":
    #         # omssa or msfragger
    #         engine = "omssa"
    #     return engine

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
                print()
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
