import bz2
import csv
import os
from decimal import Decimal
from pathlib import Path

import uparma
from chemical_composition import ChemicalComposition
from loguru import logger
from peptide_mapper.mapper import UPeptideMapper
from unimod_mapper.unimod_mapper import UnimodMapper

from unify_idents import UnifiedRow


class BaseParser:
    def __init__(self, input_file, params=None):
        pass

    @classmethod
    def file_matches_parser(self, file):
        """Check if file is compatible with parser.

        Args:
            file (str): path to file

        Returns:
            bool: Wether or not specified file can be converted by this parser.
        """
        # needs to return False to dont be selected as engine parser during `get_parsers`
        return False

    def calc_mz(self, mass, charge):
        PROTON = 1.00727646677
        return (float(mass) + (int(charge) * PROTON)) / int(charge)

    def get_column_names(self):
        headers = self.param_mapper.get_default_params(style=self.style)[
            "header_translations"
        ]["translated_value"]
        return headers

    def read_rt_lookup_file(self, scan_rt_lookup_path):
        lookup = {}
        if scan_rt_lookup_path is None:
            scan_rt_lookup_path = "\?"
        if Path(scan_rt_lookup_path).exists():
            with bz2.open(scan_rt_lookup_path, "rt") as fin:
                reader = csv.DictReader(fin)
                for line in reader:
                    file = Path(line["File"])
                    file = str(file.stem).rstrip(
                        "".join(file.suffixes)
                    )  # remove all suffixes, eg. idx.gz
                    lookup.setdefault(
                        file, {"scan2rt": {}, "rt2scan": {}, "scan2mz": {}}
                    )
                    scan, rt, mz = (
                        line["Spectrum ID"],
                        line["RT"],
                        line["Precursor mz"],
                    )
                    lookup[file]["scan2rt"][int(scan)] = float(rt)
                    lookup[file]["rt2scan"][float(rt)] = int(scan)
                    lookup[file]["scan2mz"][int(scan)] = float(mz)
        else:
            logger.warning("No scan_rt lookup file supplied")
        return lookup


class __QuantBaseParser(BaseParser):
    def __init__(self, input_file, params=None):
        super().__init__(input_file, params)
        if params is None:
            self.params = {}
        else:
            self.params = params
        self.param_mapper = uparma.UParma()
        self.cc = ChemicalComposition()
        self.scan_rt_path = self.params.get("rt_pickle_name", None)
        self.scan_rt_lookup = self.read_rt_lookup_file(self.scan_rt_path)
        self.required_headers = set(
            [
                "file_name",
                "spectrum_id",
                "trivial_name",
                "chemical_composition",
                "precursor_spectrum_id",
                "retention_time",
                "charge",
                "quant_run_id",
                "quant_value",
                "quant_score",
                "quant_group",
                "processing_level",
                "delta_mz",
                "label",
                "condition",
                "ident_reference",
                "fwhm",
                "s2i",
                "p2t",
                "coalescence",
            ]
        )

    def check_required_headers(self, row):
        return self.required_headers.issubset(row.keys())

    def __iter__(self):
        return self


class __IdentBaseParser(BaseParser):
    def __init__(self, input_file, params=None):
        super().__init__(input_file, params)
        if params is None:
            self.params = {}
        else:
            self.params = params

        self.param_mapper = uparma.UParma()
        self.mod_mapper = UnimodMapper()
        self.cc = ChemicalComposition()

        if self.params.get("modifications", None) is not None:
            self.params["mods"] = self.mod_mapper.map_mods(
                mod_list=self.params["modifications"]
            )
        else:
            self.params["mods"] = self.mod_mapper.map_mods(mod_list=[])

        self.scan_rt_path = self.params.get("rt_pickle_name", None)
        self.scan_rt_lookup = self.read_rt_lookup_file(self.scan_rt_path)
        self.create_mod_dicts()

        self.cols_to_remove = []
        self.cols_to_add = []

        no_decimals = 4
        self.mass_format_string = "{{0:3.{0}f}}".format(no_decimals)

        self.PROTON = 1.00727646677

    def __iter__(self):
        return self

    def general_fixes(self, row):
        """Apply fixed applicable to all engine parsers.

        Args:
            row (dict): dict containing psm based data from engine file

        Returns:
            dict: row with applied fixes
        """
        if row.get("Raw data location") is None or row["Raw data location"] == "":
            row["Raw data location"] = row["Spectrum Title"].split(".")[0]
        if ".mgf" in row["Raw data location"]:
            row["Raw data locations"] = row["Raw data location"].replace(
                ".mgf", ".mzML"
            )
        basename = os.path.basename(row["Raw data location"]).split(".")[0]
        if basename in self.scan_rt_lookup:
            row["Retention Time (s)"] = float(
                self.scan_rt_lookup[basename]["scan2rt"][int(row["Spectrum ID"])]
            )
            row["Exp m/z"] = self.scan_rt_lookup[basename]["scan2mz"][
                int(row["Spectrum ID"])
            ]
        else:
            row["Retention Time (s)"] = ""
            row["Exp m/z"] = ""
        row["Sequence"] = row["Sequence"].upper()

        return row

    def create_mod_dicts(self):
        """Create dict containing meta information about static and variable mods."""
        self.fixed_mods = {}
        self.opt_mods = {}
        self.mod_dict = {}
        self.n_term_replacement = {
            "Ammonia-loss": None,
            "Trimethyl": None,
            "Gly->Val": None,
        }
        self.mod_dict = {}
        self.n_term_replacement = {}
        # self.opt_mods
        for mod_type in ["fix", "opt"]:
            for modification in self.params["mods"][mod_type]:
                aa = modification["aa"]
                pos = modification["position"]
                name = modification["name"]
                if name not in self.mod_dict.keys():
                    self.mod_dict[name] = {
                        "mass": modification["mass"],
                        "aa": set(),
                        "position": set(),
                    }
                self.mod_dict[name]["aa"].add(aa)

                self.mod_dict[name]["aa"].add(pos)
                self.mod_dict[name]["position"].add(pos)

                if "N-term" in pos:
                    self.n_term_replacement[name] = aa
                if mod_type == "fix":
                    self.fixed_mods[aa] = name
                    if aa == "C" and name == "Carbamidomethyl":
                        cam = True
                        self.mod_dict["Carbamidomethyl"]["aa"].add("U")
                        self.fixed_mods["U"] = "Carbamidomethyl"
                if mod_type == "opt":
                    self.opt_mods[aa] = name

    def map_mod_names(self, row):
        """Map massshifts to unimod names.

        Args:
            row (dict): dict containing psm based data from engine file

        Returns:
            str: Description
        """
        # 0 based indexing, is corrected in this method,
        mods = []
        for mod in row["Modifications"]:
            mass, pos = mod.split(":")
            potential_names = self.mod_mapper.appMass2name_list(
                round(float(mass), 4), decimal_places=4
            )
            for name in potential_names:
                if name in self.mod_dict:
                    if row["Sequence"][int(pos)] in self.mod_dict[name]["aa"]:
                        pos = int(pos) + 1
                        mods.append(f"{name}:{pos}")  # minus
                    elif "Prot-N-term" in self.mod_dict[name]["position"]:
                        # n-term mod
                        mods.append(f"{name}:0")
        return ";".join(mods)
