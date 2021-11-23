"""Engine parser."""
import pandas as pd
import regex as re

from unify_idents.engine_parsers.base_parser import IdentBaseParser


class MSAmanda_2_Parser(IdentBaseParser):
    """File parser for MS Amanda 2."""

    def __init__(self, *args, **kwargs):
        """Initialize parser.

        Reads in data file and provides mappings.
        """
        super().__init__(*args, **kwargs)
        self.style = "msamanda_style_1"

        self.df = pd.read_csv(self.input_file, delimiter="\t", skiprows=1)
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
        # It is a csv file even though it is technically tab-delimited
        is_csv = file.as_posix().endswith(".csv")
        with open(file.as_posix()) as f:
            head = "".join([next(f) for _ in range(1)])
        matches_version = "#version: 2." in head
        return is_csv and matches_version

    def _map_mod_translation(self, row):
        """Replace single mod string.

        Args:
            row (str): unprocessed modification string

        Returns:
            mod_str (str): formatted modification string
        """
        mod_str = ""
        if row == "" or row == [""]:
            return mod_str
        for mod in row:
            mod_name = re.search(r"\(([^|]+)", mod).group(1)
            pos = mod.split("(")[0]
            if "N-TERM" in pos.upper():
                pos = 0
            else:
                pos = int(re.search(r"\d+", mod).group(0))
            mod_str += f"{mod_name}:{pos};"
        return mod_str

    def translate_mods(self):
        """
        Replace internal modification nomenclature with formatted modification strings.

        Returns:
            (pd.Series): column with formatted mod strings
        """
        mod_split_col = self.df["Modifications"].fillna("").str.split(";")
        mods_translated = mod_split_col.apply(self._map_mod_translation)

        return mods_translated.str.rstrip(";")

    def unify(self):
        """
        Primary method to read and unify engine output.

        Returns:
            self.df (pd.DataFrame): unified dataframe
        """
        self.df["Search Engine"] = "msamanda_2_0_0_17442"
        self.df["Raw data location"] = self.params["Raw data location"]
        self.df["Spectrum ID"] = self.df["Spectrum Title"].str.split(".").str[1]
        self.df["Modifications"] = self.translate_mods()
        self.process_unify_style()

        return self.df
