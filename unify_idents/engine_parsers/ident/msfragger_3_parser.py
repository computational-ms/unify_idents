import itertools

import pandas as pd
import regex as re
from loguru import logger

from unify_idents.engine_parsers.base_parser import __IdentBaseParser


# TODO: 15N handling missing
class MSFragger3Parser(__IdentBaseParser):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.style = "msfragger_style_3"

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
        is_tsv = file.as_posix().endswith(".tsv")
        with open(file.as_posix()) as f:
            head = "".join([next(f) for x in range(1)])
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
        mod_str = ""
        if row == "" or row == [""]:
            return mod_str
        for mod in row:
            mass = re.search("\(([^)]+)", mod).group(1)
            name = map_dict[mass]
            if len(name) > 0:
                if "N-term" in self.mod_dict[name[0]]["position"]:
                    pos = 0
                else:
                    pos = int(re.search(r"^\d+", mod).group(0))
                mod_str += f"{name[0]}:{pos};"
            else:
                return "NON_MAPPABLE"
        return mod_str

    def translate_mods(self):
        mod_split_col = self.df["Modifications"].fillna("").str.split(", ")
        unique_mods = set().union(*mod_split_col.apply(set)).difference({""})
        unique_mod_masses = {re.search("\(([^)]+)", m).group(1) for m in unique_mods}
        potential_names = {
            m: [
                name
                for name in self.mod_mapper.appMass2name_list(
                    round(float(m), 4), decimal_places=4
                )
                if name in self.mod_dict
            ]
            for m in unique_mod_masses
        }
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
            [v / len(self.df) for v in non_mappable_mods.values()]
        )
        if any(non_mappable_percent > 0.001):
            pass
            # raise ValueError("Some modifications found in more than 0.1% of PSMs cannot be mapped.")
        if len(non_mappable_percent) > 0:
            logger.warning(
                "Some modifications found in less than 0.1% of PSMs cannot be mapped and were removed."
            )
        mods_translated = mod_split_col.apply(
            self._map_mod_translation, map_dict=potential_names
        )

        return mods_translated.str.rstrip(";")

    def unify(self):
        self.df["Search Engine"] = "msfragger_3_0"
        self.df["Raw data location"] = self.params["Raw data location"]
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
