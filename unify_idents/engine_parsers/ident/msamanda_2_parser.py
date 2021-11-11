import pandas as pd
import regex as re

from unify_idents.engine_parsers.base_parser import __IdentBaseParser


class MSAmanda_2_Parser(__IdentBaseParser):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.style = "msamanda_style_1"

        self.df = pd.read_csv(self.input_file, delimiter="\t", skiprows=1)
        self.df.dropna(axis=1, how="all", inplace=True)

        cols_to_remove = [
            "proteinacc_start_stop_pre_post_;",
            "Filename",
            # "Rank",
        ]
        self.mapping_dict = {
            v: k
            for k, v in self.param_mapper.get_default_params(style=self.style)[
                "header_translations"
            ]["translated_value"].items()
            if k not in cols_to_remove
        }
        self.df.rename(columns=self.mapping_dict, inplace=True)
        self.df.columns = self.df.columns.str.lstrip(" ")
        if not "Modifications" in self.df.columns:
            self.df["Modifications"] = ""
        self.reference_dict.update({k: None for k in self.mapping_dict.values()})

    @classmethod
    def check_parser_compatibility(cls, file):
        is_csv = file.as_posix().endswith(".csv")
        with open(file.as_posix()) as f:
            head = "".join([next(f) for x in range(1)])
        matches_version = "#version: 2.0.0.17442" in head
        return is_csv and matches_version

    def _map_mod_translation(self, row):
        mod_str = ""
        if row == "" or row == [""]:
            return mod_str
        for mod in row:
            mod_name = re.search("\(([^|]+)", mod).group(1)
            pos = mod.split("(")[0]
            if "N-TERM" in pos.upper():
                pos = 0
            else:
                pos = int(re.search(r"\d+", mod).group(0))
            mod_str += f"{mod_name}:{pos};"
        return mod_str

    def translate_mods(self):
        mod_split_col = self.df["Modifications"].fillna("").str.split(";")
        mods_translated = mod_split_col.apply(self._map_mod_translation)

        return mods_translated.str.rstrip(";")

    def unify(self):
        self.df["Search Engine"] = "msamanda_2_0_0_17442"
        self.df["Raw data location"] = self.params["Raw data location"]
        self.df["Spectrum ID"] = self.df["Spectrum Title"].str.split(".").str[1]
        self.df["Modifications"] = self.translate_mods()
        self.process_unify_style()

        return self.df
