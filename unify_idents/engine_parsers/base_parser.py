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

        self.map_mods(self.params["modifications"])

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
        if ".mgf" in row["Raw data location"]:
            row["Raw data locations"] = row["Raw data location"].replace(
                ".mgf", ".mzML"
            )
        basename = os.path.basename(row["Raw data location"]).split(".")[0]
        row["Retention Time (s)"] = float(
            self.scan_rt_lookup[basename]["scan2rt"][int(row["Spectrum ID"])]
        )
        row["Sequence"] = row["Sequence"].upper()

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
                pos = modification["pos"]
                name = modification["name"]
                if name not in self.mod_dict.keys():
                    self.mod_dict[name] = {
                        "mass": modification["mass"],
                        "aa": set(),
                        "pos": set(),
                    }
                self.mod_dict[name]["aa"].add(aa)

                self.mod_dict[name]["aa"].add(pos)
                self.mod_dict[name]["pos"].add(pos)

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
                    elif "Prot-N-term" in self.mod_dict[name]["pos"]:
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

    def map_mods(self, mods):
        # TODO remove logger.warning functions and replace by logger
        self.params["mods"] = {"fix": [], "opt": []}
        for ursgal_index, mod in enumerate(sorted(mods)):
            mod_params = mod.split(",")
            if len(mod_params) >= 6 or len(mod_params) <= 3:
                logger.warning(
                    """
For modifications, please use the ursgal_style:
'amino_acid,opt/fix,position,Unimod PSI-MS Name'
or
'amino_acid,opt/fix,position,name,chemical_composition'
Continue without modification {0} """.format(
                        mod
                    )
                )
                continue
            aa = mod_params[0].strip()
            mod_option = mod_params[1].strip()
            pos = mod_params[2].strip()
            unimod = False
            unimod_id = None

            if len(mod_params) == 4:
                try:
                    unimod_id = int(mod_params[3].strip())
                    unimod_name = self.mod_mapper.id2name(unimod_id)
                    mass = self.mod_mapper.id2mass(unimod_id)
                    composition = self.mod_mapper.id2composition(unimod_id)
                    if unimod_name is None:
                        logger.warning(
                            """
'{1}' is not a Unimod modification
please change it to a valid Unimod Accession # or PSI-MS Unimod Name
or add the chemical composition hill notation (including 1)
e.g.: H-1N1O2
ursgal_style: 'amino_acid,opt/fix,position,name,chemical_composition'
Continue without modification {0} """.format(
                                mod, unimod_id
                            )
                        )
                        continue
                    unimod = True
                    name = unimod_name
                except:
                    unimod_name = mod_params[3].strip()
                    unimod_id = self.mod_mapper.name2id(unimod_name)
                    mass = self.mod_mapper.name2mass(unimod_name)
                    composition = self.mod_mapper.name2composition(unimod_name)
                    if unimod_id is None:
                        logger.warning(
                            """
'{1}' is not a Unimod modification
please change it to a valid PSI-MS Unimod Name or Unimod Accession #
or add the chemical composition hill notation (including 1)
e.g.: H-1N1O2
ursgal_style: 'amino_acid,opt/fix,position,name,chemical_composition'
Continue without modification {0} """.format(
                                mod, unimod_name
                            )
                        )
                        continue
                    unimod = True
                    name = unimod_name

            elif len(mod_params) == 5:
                name = mod_params[3].strip()
                chemical_formula = mod_params[4].strip()
                chemical_composition = ursgal.ChemicalComposition()
                chemical_composition.add_chemical_formula(chemical_formula)
                composition = chemical_composition
                composition_unimod_style = chemical_composition.hill_notation_unimod()
                unimod_name_list = self.mod_mapper.composition2name_list(
                    composition_unimod_style
                )
                unimod_id_list = self.mod_mapper.composition2id_list(
                    composition_unimod_style
                )
                mass = self.mod_mapper.composition2mass(composition_unimod_style)
                for i, unimod_name in enumerate(unimod_name_list):
                    if unimod_name == name:
                        unimod_id = unimod_id_list[i]
                        unimod = True
                        break
                if unimod == False and unimod_name_list != []:
                    logger.warning(
                        """
'{0}' is not a Unimod modification
but the chemical composition you specified is included in Unimod.
Please use one of the Unimod names:
{1}
Continue without modification {2} """.format(
                            name, unimod_name_list, mod
                        )
                    )
                    continue
                if unimod == False and unimod_name_list == []:
                    logger.warning(
                        """
'{0}' is not a Unimod modification
trying to continue with the chemical composition you specified
This is not working with OMSSA so far""".format(
                            mod,
                        )
                    )
                    mass = chemical_composition._mass()
                    # write new userdefined modifications Xml in unimod style

            mod_dict = {
                "_id": ursgal_index,
                "aa": aa,
                "mass": mass,
                "pos": pos,
                "name": name,
                "composition": composition,
                "org": mod,
                "id": unimod_id,
                "unimod": unimod,
            }
            if mod_dict["unimod"] == False:
                self.mod_mapper.writeXML(mod_dict)

            self.params["mods"][mod_option].append(mod_dict)
        return self.params["mods"]
