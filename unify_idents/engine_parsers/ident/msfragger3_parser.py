import csv
import itertools
import re
from decimal import ROUND_UP, Decimal, getcontext
from pathlib import Path

from loguru import logger

from unify_idents import UnifiedRow
from unify_idents.engine_parsers.base_parser import __IdentBaseParser

"""
1. Raw data location
2. Mass and m/z calculation
3. Remapping proteins
4. 

"""


class MSFragger3Parser(__IdentBaseParser):

    """Initialize MSFragger parser.

    Args:
        input_file (str): path to file to unify
        params (dict, optional): parser specific parameters
    """

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
        self.mod_pattern = re.compile(r""":(?P<pos>[0-9]*$)""")
        self.mod_pattern_msfragger_term = re.compile(
            r""".-term\((?P<mass>[0-9]*\.[0-9]*)\)"""
        )
        self.mod_pattern_msfragger = re.compile(
            r"""(?P<pos>[0-9]*)(?P<aa>[A-Z])\((?P<mass>[0-9]*\.[0-9]*)\)"""
        )
        self.mass_to_mod_combo = self.prepare_mass_to_mod()

    def __next__(self):
        """Return next unified line from reader.

        Returns:
            UnifiedRow: unified PSM level data.
        """
        line = next(self.reader)
        line = self._unify_row(line)
        return line

    @classmethod
    def file_matches_parser(cls, file):
        """Check if file is compatible with parser.

        Args:
            file (str): path to file

        Returns:
            bool: Wether or not specified file can be converted by this parser.
        """
        column_names = [
            "scannum",
            "peptide",
            "charge",
            "peptide_prev_aa",
            "peptide_next_aa",
            "protein",
            "modification_info",
            "retention_time",
            "precursor_neutral_mass",
            "calc_neutral_pep_mass",
            "hit_rank",
            "massdiff",
            "num_matched_ions",
            "tot_num_ions",
            "hyperscore",
            "nextscore",
            "num_tol_term",
            "num_missed_cleavages",
            "expectscore",
            "best_locs",
            "score_without_delta_mass",
            "best_score_with_delta_mass",
            "second_best_score_with_delta_mass",
            "delta_score",
        ]
        with open(file) as fin:
            headers = fin.readline()
            if set(headers.split()) == set(column_names):
                ret_val = True
            else:
                ret_val = False
        return ret_val

    def _unify_row(self, row):
        """Convert row to unified format.

        Args:
            row (dict): dict containing psm based ident information.

        Returns:
            UnifiedRow: converted row
        """
        new_row = {}
        for unify_name, omssa_name in self.column_mapping.items():
            new_row[unify_name] = row[omssa_name]
        for col in self.cols_to_remove:
            del new_row[col]
        for col in self.cols_to_add:
            if col not in new_row:
                new_row[col] = ""
        new_row["Search Engine"] = "msfragger_3_0"
        new_row["Spectrum Title"] = "{file}.{specid}.{specid}.{charge}".format(
            file=str(self.params["Raw data location"]).split(".")[0],
            specid=new_row["Spectrum ID"],
            charge=new_row["Charge"],
        )
        new_row["Raw data location"] = self.params["Raw data location"]
        new_row["Exp m/z"] = self.calc_mz(
            new_row["MSFragger:Precursor neutral mass (Da)"], new_row["Charge"]
        )
        new_row["Calc m/z"] = self.calc_mz(
            new_row["MSFragger:Neutral mass of peptide"], new_row["Charge"]
        )

        # TODO
        # think of all the labile mode, glycan and 15N stuff ...
        modstring = self.format_mods(new_row)

        new_row["Modifications"] = modstring

        new_row = self.general_fixes(new_row)
        return UnifiedRow(**new_row)

    def prepare_mass_to_mod(self):
        """Map massshifts to unimod names.

        Args:
            row (dict): dict containing psm based data from engine file

        Returns:
            str: Description
        """
        # rounding using Decimal (stdlib)
        getcontext().prec = 8
        getcontext().rounding = ROUND_UP
        # mod_dict_list = params['mods']['opt'] + params['mods']['fix']

        # 1.2 Add Mod 15N combinations to self.mod_dict
        # TODO get rid of use15N param
        use15N = self.params.get("15N", False)
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

        # 2. Create all combinations of all possible mods, calc masses for these "merged mods" and save a mass to mod combo lookup
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
        """Convert mods to unified modstring.

        Args:
            new_row (dict): unified row

        Returns:
            str: unified modstring
        """
        final_mods = []
        for single_mod in row["Modifications"].split(", "):
            if single_mod.strip() == "":
                continue
            else:
                match = self.mod_pattern_msfragger.search(single_mod)
                pos = int(match.group("pos"))
                mass = float(match.group("mass"))
                aa = match.group("aa")
                mass_string = msfragger_mass = self.mass_format_string.format(
                    Decimal(mass)
                )

            if mass_string in self.mass_to_mod_combo.keys():
                explainable_combos = []
                for combo in self.mass_to_mod_combo[mass_string]:
                    for new_name in combo:
                        meta_mod_info = self.mod_dict[new_name]
                        real_pos = pos
                        if (
                            "Prot-N-term" in self.mod_dict[new_name]["position"]
                            or "N-term" in self.mod_dict[new_name]["position"]
                        ):
                            real_pos = 0  # should be 1 and is set to 0
                        final_mods.append(f"{new_name}:{real_pos}")
            else:
                # not a merged mod
                names = self.mod_mapper.appMass2name_list(mass, decimal_places=6)
                for new_name in names:
                    if new_name in self.mod_dict:
                        real_pos = pos
                        if (
                            "Prot-N-term" in self.mod_dict[new_name]["position"]
                            or "N-term" in self.mod_dict[new_name]["position"]
                        ):
                            real_pos = 0  # should be 1 and is set to 0
                        final_mods.append(f"{new_name}:{pos}")
                        break
        return ";".join(final_mods)
