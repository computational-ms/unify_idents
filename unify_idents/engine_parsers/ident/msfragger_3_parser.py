"""Engine parser."""
import itertools

import pandas as pd
import regex as re
from loguru import logger

from unify_idents.engine_parsers.base_parser import IdentBaseParser


class MSFragger_3_Parser(IdentBaseParser):
    """File parser for MSFragger 3."""

    def __init__(self, *args, **kwargs):
        """Initialize parser.

        Reads in data file and provides mappings.
        """
        super().__init__(*args, **kwargs)
        self.style = "msfragger_style_3"
        # 15N handling missing for now
        if self.params.get("label", "") == "15N":
            raise NotImplementedError

        self.df = pd.read_csv(self.input_file, delimiter="\t")
        self.df.dropna(axis=1, how="all", inplace=True)

        self.mapping_dict = {
            v: k
            for k, v in self.param_mapper.get_default_params(style=self.style)[
                "header_translations"
            ]["translated_value"].items()
        }
        self.df.rename(columns=self.mapping_dict, inplace=True)
        self.df.columns = self.df.columns.str.lstrip(" ")
        if not "Modifications" in self.df.columns:
            self.df["Modifications"] = ""
        self.reference_dict.update({k: None for k in self.mapping_dict.values()})

    @classmethod
    def check_parser_compatibility(cls, file):
        """Assert compatibility between file and parser.

        Args:
            file (str): path to input file

        Returns:
            bool: True if parser and file are compatible

        """
        is_tsv = file.as_posix().endswith(".tsv")
        with open(file.as_posix()) as f:
            head = "".join([next(f) for _ in range(1)])
        head = set(head.rstrip("\n").split("\t"))
        ref_columns = {
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
        }
        columns_match = len(ref_columns.difference(head)) == 0
        return is_tsv and columns_match

    def _map_mod_translation(self, row, map_dict):
        """Replace single mod string.

        Args:
            row (str): unprocessed modification string
            map_dict (dict): mod mapping dict

        Returns:
            mod_str (str): formatted modification string
        """
        mod_str = ""
        if row == "" or row == [""]:
            return mod_str
        for mod in row:
            mass = re.search(r"\(([^)]+)", mod).group(1)
            name = map_dict[mass]
            if len(name) > 0:
                for m in name:
                    if any(["N-term" in p for p in self.mod_dict[m]["position"]]):
                        pos = 0
                    else:
                        pos = int(re.search(r"^\d+", mod).group(0))
                    mod_str += f"{m}:{pos};"
            else:
                return "NON_MAPPABLE"
        return mod_str

    def translate_mods(self):
        """
        Replace internal modification nomenclature with formatted modification strings.

        Returns:
            (pd.Series): column with formatted mod strings
        """
        mod_split_col = self.df["Modifications"].fillna("").str.split(", ")
        unique_mods = set().union(*mod_split_col.apply(set)).difference({""})
        unique_mod_masses = {re.search(r"\(([^)]+)", m).group(1) for m in unique_mods}
        # Map single mods
        potential_names = {
            m: [
                name
                for name in self.mod_mapper.mass_to_names(round(float(m), 4), decimals=4)
                if name in self.mod_dict
            ]
            for m in unique_mod_masses
        }
        # Map multiple mods
        for n in [2, 3]:
            for unmapped_mass in {k: v for k, v in potential_names.items() if v == []}:
                potential_mods = [
                    name[1]
                    for name in self.mod_mapper.mass_to_combos(
                        round(float(unmapped_mass), 4), n=n, decimals=4
                    )
                    if all(m in self.mod_dict for m in name[1])
                ]
                if len(potential_mods) == 1:
                    potential_names[unmapped_mass] = potential_mods[0]
        non_mappable_mods = {
            k: len(
                [
                    m
                    for m in list(
                        itertools.chain.from_iterable(
                            mod_split_col.apply(list).to_list()
                        )
                    )
                    if k in m
                ]
            )
            for k, v in potential_names.items()
            if v == []
        }
        non_mappable_percent = pd.Series(
            [v / len(self.df) for v in non_mappable_mods.values()], dtype="float64"
        )
        if any(non_mappable_percent > 0.001):
            raise ValueError(
                "Some modifications found in more than 0.1% of PSMs cannot be mapped."
            )
        if len(non_mappable_percent) > 0:
            logger.warning(
                "Some modifications found in less than 0.1% of PSMs cannot be mapped and were removed."
            )
        mods_translated = mod_split_col.apply(
            self._map_mod_translation, map_dict=potential_names
        )

        return mods_translated.str.rstrip(";")

    def unify(self):
        """
        Primary method to read and unify engine output.

        Returns:
            self.df (pd.DataFrame): unified dataframe
        """
        self.df["Search Engine"] = "msfragger_3_0"
        self.df["Raw data location"] = str(self.params["Raw data location"])
        spec_title = (
            self.df["Raw data location"].str.split(".").str[0].str.split("/").str[-1]
        )
        self.df["Spectrum Title"] = (
            spec_title
            + "."
            + self.df["Spectrum ID"].astype(str)
            + "."
            + self.df["Charge"].astype(str)
        )
        self.df["Exp m/z"] = self._calc_mz(
            mass=self.df["MSFragger:Precursor neutral mass (Da)"],
            charge=self.df["Charge"],
        )
        self.df["Calc m/z"] = self._calc_mz(
            mass=self.df["MSFragger:Neutral mass of peptide"],
            charge=self.df["Charge"],
        )
        self.df["Modifications"] = self.translate_mods()
        self.df = self.df.loc[
            ~self.df["Modifications"].str.contains("NON_MAPPABLE", regex=False), :
        ]
        self.process_unify_style()

        return self.df
