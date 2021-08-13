import bz2
import csv
from pathlib import Path
import os

import uparma
from peptide_mapper.mapper import UPeptideMapper
from unimod_mapper.unimod_mapper import UnimodMapper
from chemical_composition import ChemicalComposition

from unify_idents import UnifiedRow
from loguru import logger
from decimal import Decimal

# get from params


class __BaseParser:
    def __init__(self, input_file, params=None):
        if params is None:
            self.params = {}
        else:
            self.params = params

        self.param_mapper = uparma.UParma()
        # self.peptide_mapper = UPeptideMapper(params["database"])
        self.mod_mapper = UnimodMapper()
        self.cc = ChemicalComposition()

        # self.map_mods(self.params["modifications"])
        self.params["mods"] = self.mod_mapper.map_mods(mod_list=self.params["modifications"])
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

    @classmethod
    def file_matches_parser(self, file):
        # needs to return False to dont be selected as engine parser during `get_parsers`
        return False

    def general_fixes(self, row):
        if row.get("Raw data location") is None or row["Raw data location"] == "":
            row["Raw data location"] = row["Spectrum Title"].split(".")[0]
        # if ".mgf" in str(row["Raw data location"]):
        #     row["Raw data locations"] = row["Raw data location"].str.replace(
        #         ".mgf", ".mzML"
        #     )
        basename = os.path.basename(row["Raw data location"]).split(".")[0]
        row["Retention Time (s)"] = float(
            self.scan_rt_lookup[basename]["scan2rt"][int(row["Spectrum ID"])]
        )
        row["Sequence"] = row["Sequence"].upper()
        row["Exp m/z"] = self.scan_rt_lookup[basename]["scan2mz"][
            int(row["Spectrum ID"])
        ]

        return row

    def check_mod_positions(self, row):
        return row

    # currently not used
    # def recalc_masses(row):
    #     self.cc.use(sequence=row["Sequence"], modifications=row["Modifications"])
    #     row["uCalc m/z"] = self.calc_mz(self.cc.mass(), int(row["Charge"]))
    #     row["uCalc mass"] = self.cc.mass()
    #     return row

    def calc_mz(self, mass, charge):
        PROTON = 1.00727646677
        return (float(mass) + (int(charge) * PROTON)) / int(charge)

    def create_mod_dicts(self):
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

    def map_peptides(self, row):
        starts = []
        ids = []
        stops = []
        pre = []
        post = []
        # we need to convert sequences to uppercase, e.g. omssa reports modified AAs in lowercase
        mapped = self.peptide_mapper.map_peptides(
            [row["Sequence"].upper()]
        )  # uses 99% of time
        for seq, data_list in mapped.items():
            for data in data_list:
                ids.append(data["id"])
                starts.append(str(data["start"]))
                stops.append(str(data["end"]))
                pre.append(str(data["pre"]))
                post.append(str(data["post"]))
        row["Protein ID"] = DELIMITER.join(ids)
        row["Sequence Pre AA"] = DELIMITER.join(pre)
        row["Sequence Post AA"] = DELIMITER.join(post)
        row["Sequence Start"] = DELIMITER.join(starts)
        row["Sequence Stop"] = DELIMITER.join(stops)
        return row

    def read_rt_lookup_file(self, scan_rt_lookup_path):
        with bz2.open(scan_rt_lookup_path, "rt") as fin:
            lookup = {}
            reader = csv.DictReader(fin)
            for line in reader:
                file = Path(line["File"])
                file = str(file.stem).rstrip(
                    "".join(file.suffixes)
                )  # remove all suffixes, eg. idx.gz
                lookup.setdefault(file, {"scan2rt": {}, "rt2scan": {}, "scan2mz": {}})
                scan, rt, mz = (
                    line["Spectrum ID"],
                    line["RT"],
                    line["Precursor mz"],
                )
                lookup[file]["scan2rt"][int(scan)] = float(rt)
                lookup[file]["rt2scan"][float(rt)] = int(scan)
                lookup[file]["scan2mz"][int(scan)] = float(mz)
        return lookup
