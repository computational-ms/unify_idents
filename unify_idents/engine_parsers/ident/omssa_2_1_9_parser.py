"""Engine parser."""

import numpy as np
import pandas as pd

from unify_idents.engine_parsers.base_parser import IdentBaseParser


class Omssa_Parser(IdentBaseParser):
    """File parser for OMSSA."""

    def __init__(self, *args, **kwargs):
        """Initialize parser.

        Reads in data file and provides mappings.
        """
        super().__init__(*args, **kwargs)
        self.style = "omssa_style_1"

        self.df = pd.read_csv(self.input_file)

        self.mapping_dict = {
            v: k
            for k, v in self.param_mapper.get_default_params(style=self.style)[
                "header_translations"
            ]["translated_value"].items()
        }
        self.df.rename(columns=self.mapping_dict, inplace=True)
        self.df.columns = self.df.columns.str.lstrip(" ")
        self.df.drop(
            columns=[
                c
                for c in self.df.columns
                if c
                not in set(self.mapping_dict.values()) | set(self.reference_dict.keys())
            ],
            inplace=True,
            errors="ignore",
        )
        self.reference_dict.update({k: None for k in self.mapping_dict.values()})

    @classmethod
    def check_parser_compatibility(cls, file):
        """Assert compatibility between file and parser.

        Args:
            file (str): path to input file

        Returns:
            bool: True if parser and file are compatible

        """
        is_csv = file.as_posix().endswith(".csv")
        with open(file.as_posix()) as f:
            try:
                head = "".join([next(f) for _ in range(1)])
            except StopIteration:
                head = ""
        head = set(head.rstrip("\n").split(","))
        ref_columns = {
            "Spectrum number",
            " Filename/id",
            " Peptide",
            " E-value",
            " Mass",
            " gi",
            " Accession",
            " Start",
            " Stop",
            " Defline",
            " Mods",
            " Charge",
            " Theo Mass",
            " P-value",
            " NIST score",
        }
        columns_match = len(ref_columns.difference(head)) == 0
        return is_csv and columns_match

    def _replace_mod_strings(self, row, mod_translations):
        """Replace single mod string.

        Args:
            row (str): unprocessed modification string
            mod_translations (dict): mod translation dict

        Returns:
            mod_str (str): formatted modification string
        """
        if row == "":
            return ""
        mods = row.split(" ,")
        new_modstring = ""
        for m in mods:
            omssa_name, pos = m.split(":")
            unimod_name = mod_translations[omssa_name]["unimod_name"]
            if pos == "1" and any(
                [
                    ("N-term" in target)
                    for target in mod_translations[omssa_name]["aa_targets"]
                ]
            ):
                pos = "0"
            new_modstring += f"{unimod_name}:{pos};"
        return new_modstring.rstrip(";")

    def translate_mods(self):
        """
        Replace internal modification nomenclature with formatted modification strings.

        Operations are performed inplace.
        """
        self.df["sequence"] = self.df["sequence"].str.upper()
        self.df["modifications"] = (
            self.df["modifications"].fillna("").str.replace(" ,", ";")
        )
        fix_mods = None
        # Map fixed mods
        fixed_mod_types = [
            d for d in self.params["modifications"] if d["type"] == "fix"
        ]
        for fm in fixed_mod_types:
            fm_strings = (
                self.df["sequence"]
                .str.split(fm["aa"])
                .apply(
                    lambda l: ";".join(
                        [
                            fm["name"] + ":" + ind
                            for ind in (
                                np.cumsum(list(map(len, l[:-1]))) + range(1, len(l))
                            ).astype(str)
                        ]
                    )
                )
            )
            if fix_mods is None:
                fix_mods = fm_strings
            else:
                fix_mods = fix_mods + ";" + fm_strings

        if fix_mods is not None:
            self.df["modifications"] = (
                self.df["modifications"].fillna("") + ";" + fix_mods
            )

    def unify(self):
        """
        Primary method to read and unify engine output.

        Returns:
            self.df (pd.DataFrame): unified dataframe
        """
        self.df["calc_mz"] = self._calc_mz(
            mass=self.df["calc_mz"], charge=self.df["charge"]
        )
        self.df["spectrum_id"] = (
            self.df["spectrum_title"].str.split(".").str[-3].astype(int)
        )
        self.df["search_engine"] = "omssa_2_1_9"
        self.translate_mods()
        self.process_unify_style()

        return self.df
