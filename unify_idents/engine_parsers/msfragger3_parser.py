from unify_idents import UnifiedRow
from unify_idents.engine_parsers.base_parser import __BaseParser
from pathlib import Path
import re
import csv
from decimal import Decimal, getcontext, ROUND_UP
import itertools

"""
1. Raw data location
2. Mass and m/z calculation
3. Remapping proteins
4. 

"""


class MSFragger3Parser(__BaseParser):
    def __init__(self, input_file, params=None):
        super().__init__(input_file, params)
        if params is None:
            params = {}
        self.params = params
        self.input_file = input_file

        try:
            self.reader = csv.DictReader(open(input_file), delimiter="\t")
        except:
            self.reader = iter([])

        self.style = "msfragger_style_3"
        self.column_mapping = self.get_column_names()

        self.cols_to_add = [
            "Raw file location",
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

        self.DICT_15N_DIFF = {
            "A": 0.997035,
            "C": 0.997035,
            "D": 0.997035,
            "E": 0.997035,
            "F": 0.997035,
            "G": 0.997035,
            "H": 2.991105,
            "I": 0.997035,
            "K": 1.99407,
            "L": 0.997035,
            "M": 0.997035,
            "N": 1.99407,
            "P": 0.997035,
            "Q": 1.99407,
            "R": 3.98814,
            "S": 0.997035,
            "T": 0.997035,
            "V": 0.997035,
            "W": 1.99407,
            "Y": 0.997035,
        }

        self.mass2mod = {str(d["mass"]): d["name"] for d in self.params["mods"]["opt"]}
        self.mass2mod.update(
            {str(d["mass"]): d["name"] for d in self.params["mods"]["fix"]}
        )
        self.mod_pattern_msfragger_term = re.compile(
            r""".-term\((?P<mass>[0-9]*\.[0-9]*)\)"""
        )
        self.mod_pattern_msfragger = re.compile(
            r"""(?P<pos>[0-9]*)(?P<aa>[A-Z])\((?P<mass>[0-9]*\.[0-9]*)\)"""
        )
        no_decimals = 5
        self.mass_format_string = "{{0:3.{0}f}}".format(no_decimals)
        self.mass_to_mod_combo = self.prepare_mass_to_mod()

    def __next__(self):
        line = next(self.reader)
        line = self._unify_row(line)
        return line

    def file_matches_parser(self):
        with open(self.input_file) as fin:
            headers = fin.readline()
            if set(headers.split()) == set(self.column_mapping.values()):
                ret_val = True
            else:
                ret_val = False
        return ret_val

    def _unify_row(self, row):
        new_row = {}
        for unify_name, omssa_name in self.column_mapping.items():
            new_row[unify_name] = row[omssa_name]
        for col in self.cols_to_remove:
            del new_row[col]
        for col in self.cols_to_add:
            if col not in new_row:
                # print(f"ADD {col}")
                new_row[col] = ""

        new_row["Search Engine"] = "msfragger_3_0"
        new_row["Spectrum Title"] = "{file}.{specid}.{specid}.{charge}".format(
            file=self.params["Raw file location"],
            specid=new_row["Spectrum ID"],
            charge=new_row["Charge"],
        )
        new_row["Raw file location"] = self.params["Raw file location"]
        new_row["Exp m/z"] = self.calc_mz(
            new_row["MSFragger:Precursor neutral mass (Da)"], new_row["Charge"]
        )
        new_row["Calc m/z"] = self.calc_mz(
            new_row["MSFragger:Neutral mass of peptide"], new_row["Charge"]
        )

        # TODO
        # think of all the labile mode, glycan and 15N stuff ...
        modstring = self.format_mods(new_row)
        # if modstring != "":
        #     breakpoint()

        new_row["Modifications"] = modstring

        new_row = self.general_fixes(new_row)
        return UnifiedRow(**new_row)

    def get_column_names(self):
        headers = self.param_mapper.get_default_params(style=self.style)[
            "header_translations"
        ]["translated_value"]
        return headers

    def prepare_mass_to_mod(self):
        fixed_mods = {}
        opt_mods = {}
        self.mod_dict = {}

        n_term_replacement = {
            "Ammonia-loss": None,
            "Trimethyl": None,
            "Gly->Val": None,
        }

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
                if "N-term" in pos:
                    n_term_replacement[name] = aa
                if mod_type == "fix":
                    fixed_mods[aa] = name
                    if aa == "C" and name == "Carbamidomethyl":
                        cam = True
                        # allow also Carbamidomnethyl on U, since the mod name gets changed
                        # already in upeptide_mapper
                        # According to unimod, the modification is also on Selenocystein
                        # otherwise we should change that back so that it is skipped...
                        self.mod_dict["Carbamidomethyl"]["aa"].add("U")
                        fixed_mods["U"] = "Carbamidomethyl"
                if mod_type == "opt":
                    opt_mods[aa] = name

        getcontext().prec = 8
        getcontext().rounding = ROUND_UP
        # mod_dict_list = params['mods']['opt'] + params['mods']['fix']
        use15N = True
        if use15N:
            aminoacids_2_check = set()
            for modname in self.mod_dict.keys():
                aminoacids_2_check |= self.mod_dict[modname]["aa"]
            additional_15N_modifications = []
            for aminoacid, N15_Diff in self.DICT_15N_DIFF.items():
                if aminoacid not in aminoacids_2_check:
                    continue
                if "_15N_{0}".format(aminoacid) in self.mod_dict.keys():
                    print(
                        """
                        Error in unify_csv
                        New mod_name already present in mod_dict'
                        This should never happen"""
                    )
                    sys.exit(1)
                self.mod_dict["_15N_{0}".format(aminoacid)] = {
                    "mass": N15_Diff,
                    "aa": set([aminoacid]),
                    "pos": set(["any"]),
                }

        mod_names = []
        for mod in sorted(list(self.mod_dict.keys())):
            mod_names.extend(itertools.repeat(mod, len(self.mod_dict[mod]["aa"])))
        mass_to_mod_combo = {}
        # we cover all combocs of mods
        for iter_length in range(2, len(mod_names) + 1):
            for name_combo in itertools.combinations(mod_names, iter_length):
                mass = 0
                for name in name_combo:
                    mass += Decimal(self.mod_dict[name]["mass"])
                rounded_mass = self.mass_format_string.format(mass)
                if rounded_mass not in mass_to_mod_combo.keys():
                    mass_to_mod_combo[rounded_mass] = set()
                mass_to_mod_combo[rounded_mass].add(name_combo)
        return mass_to_mod_combo

    def format_mods(self, row):
        ms_fragger_reformatted_mods = []
        if row["Modifications"] == "M":
            # M stand for Modifications here, not Methionine
            row["Modifications"] = ""
        else:
            mod_string = row["Modifications"]
            if mod_string == "":
                mod_list = []
            elif "|" in mod_string:
                mod_list = []
                for single_mod in mod_string.split("|"):
                    if single_mod in ["M", ""]:
                        continue
                    msfragger_pos, raw_msfragger_mass = single_mod.split("$")
                    msfragger_mass = self.mass_format_string.format(
                        # mass rounded as defined above
                        Decimal(raw_msfragger_mass)
                    )
                    msfragger_pos = int(msfragger_pos)
                    mod_list.append((msfragger_mass, raw_msfragger_mass, msfragger_pos))
            else:
                mod_list = []
                for single_mod in mod_string.split(", "):
                    if single_mod.startswith("N-term"):
                        msfragger_pos = 0
                        match = self.mod_pattern_msfragger_term.search(single_mod)
                    elif single_mod.startswith("C-term"):
                        msfragger_pos = len(row["Sequence"])
                    else:
                        match = self.mod_pattern_msfragger.search(single_mod)
                        msfragger_pos = int(match.group("pos"))
                        if msfragger_pos != 0:
                            msfragger_pos -= 1
                    raw_msfragger_mass = float(match.group("mass"))
                    msfragger_mass = self.mass_format_string.format(
                        # mass rounded as defined above
                        Decimal(raw_msfragger_mass)
                    )
                    mod_list.append((msfragger_mass, raw_msfragger_mass, msfragger_pos))
            for single_mod in mod_list:
                (
                    msfragger_mass,
                    raw_msfragger_mass,
                    msfragger_pos,
                ) = single_mod
                if msfragger_mass in self.mass_to_mod_combo.keys():
                    explainable_combos = []
                    for combo in self.mass_to_mod_combo[msfragger_mass]:
                        combo_explainable = set([True])
                        tmp_mods = []
                        for new_name in combo:
                            meta_mod_info = self.mod_dict[new_name]
                            single_mod_check = set([True])
                            # check aa
                            if (
                                "*" not in meta_mod_info["aa"]
                                and row["Sequence"][msfragger_pos]
                                not in meta_mod_info["aa"]
                            ):
                                single_mod_check.add(False)
                            # check pos
                            if "any" not in meta_mod_info["pos"]:
                                pos_to_check = set()
                                if (
                                    "Prot-N-term" in meta_mod_info["pos"]
                                    or "N-term" in meta_mod_info["pos"]
                                ):
                                    pos_to_check.add(0)
                                elif (
                                    "Prot-C-term" in meta_mod_info["pos"]
                                    or "C-term" in meta_mod_info["pos"]
                                ):
                                    pos_to_check.add(int(len(row["Sequence"])) - 1)
                                else:
                                    pass
                                if pos_to_check != set():
                                    if msfragger_pos not in pos_to_check:
                                        single_mod_check.add(False)

                            if all(single_mod_check):
                                pos_in_peptide_for_format_str = msfragger_pos + 1
                                # we keep mass here so that the
                                # correct name is added later in already
                                # existing code
                                tmp_mods.append(
                                    "{0}:{1}".format(
                                        meta_mod_info["mass"],
                                        pos_in_peptide_for_format_str,
                                    )
                                )
                            else:
                                combo_explainable.add(False)
                        if all(combo_explainable):
                            explainable_combos.append(tmp_mods)
                    if len(explainable_combos) > 1:
                        print(
                            """
                            [ WARNING ] Multiple modification combinations possible
                            [ WARNING ] to explain reported modification mass
                            [ WARNING ] The following combination was chosen to continue:
                            [ WARNING ] {0}
                            """.format(
                                sorted(explainable_combos)[0],
                            )
                        )
                    elif len(explainable_combos) == 1:
                        ms_fragger_reformatted_mods += sorted(explainable_combos)[0]
                    else:
                        # no combos explainable
                        ms_fragger_reformatted_mods.append(
                            "{0}:{1}".format(raw_msfragger_mass, msfragger_pos + 1)
                        )
                else:
                    # MS Frager starts counting at zero
                    ms_fragger_reformatted_mods.append(
                        "{0}:{1}".format(raw_msfragger_mass, msfragger_pos + 1)
                    )
        # mods = ";".join(ms_fragger_reformatted_mods)
        formatted_mods = []
        for i, aa in enumerate(row["Sequence"]):
            i += 1
            for mod in self.params["mods"]["fix"]:
                if aa in mod["aa"]:
                    formatted_mods.append(f"{mod['name']}:{i}")
        for mod in ms_fragger_reformatted_mods:
            mass, pos = mod.split(":")
            mass, pos = float(mass), int(pos)
            names = self.mod_mapper.appMass2name_list(mass, decimal_places=6)
            for n in names:
                if n in self.mod_dict:
                    if row["Sequence"][pos - 1] in self.mod_dict[n]["aa"]:
                        formatted_mods.append(f"{n}:{pos}")
        formatted_mods = ";".join(formatted_mods)
        return formatted_mods
