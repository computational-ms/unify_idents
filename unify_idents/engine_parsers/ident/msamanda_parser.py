#!/usr/bin/env python
import csv

from unify_idents import UnifiedRow
from unify_idents.engine_parsers.base_parser import __IdentBaseParser


class MSamandaParser(__IdentBaseParser):

    """Engine parser to unify MSAmanda results."""

    def __init__(self, input_file, params=None):
        """Initialize MSAmanda parser.

        Args:
            input_file (str): path to file to unify
            params (dict, optional): parser specific parameters
        """
        super().__init__(input_file, params)
        if params is None:
            params = {}
        self.params = params
        self.input_file = input_file

        try:
            result_file = open(input_file, "r")
            self.reader = csv.DictReader(
                (row for row in result_file if not row.startswith("#")), delimiter="\t"
            )
        except:
            self.reader = None

        self.style = "msamanda_style_1"
        self.column_mapping = self.get_column_names()
        self.cols_to_remove = [
            "proteinacc_start_stop_pre_post_;",
            "Filename",
            # "Rank",
        ]

        self.cols_to_add = [
            "uCalc m/z",
            "uCalc Mass",
            "Calc m/z",
            "Exp m/z",
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

    @classmethod
    def file_matches_parser(cls, file):
        """Check if file is compatible with parser.

        Args:
            file (str): path to file

        Returns:
            bool: Wether or not specified file can be converted by this parser.
        """
        msamanda_version = "#version: 2.0.0.17442"
        with open(file) as fh:
            reader = csv.DictReader(fh)
            if reader.fieldnames[0] == msamanda_version:
                ret_val = True
            else:
                ret_val = False
        return ret_val

    def __iter__(self):
        return self

    def __next__(self):
        n = next(self.reader)
        u = self._unify_row(n)
        return u

    def _unify_row(self, row):
        """Convert row to unified format.

        Args:
            row (dict): dict containing psm based ident information.

        Returns:
            UnifiedRow: converted row
        """
        new_row = {}
        for unify_name, engine_name in self.column_mapping.items():
            new_row[unify_name] = row[engine_name]
        for col in self.cols_to_remove:
            del new_row[col]
        for col in self.cols_to_add:
            new_row[col] = ""
        new_row["Spectrum ID"] = int(new_row["Spectrum Title"].split(".")[1])
        new_row["Search Engine"] = "msamanda_2_0_0_17442"

        modstring = self.create_mod_string(new_row)
        new_row["Modifications"] = modstring
        new_row = self.general_fixes(new_row)

        return UnifiedRow(**new_row)

    def create_mod_string(self, new_row):
        """Convert mods to unified modstring.

        Args:
            new_row (dict): unified row

        Returns:
            str: unified modstring
        """
        mod_input = new_row["Modifications"]
        translated_mods = []
        # N-Term(Acetyl|42.010565|fixed);M1(Oxidation|15.994915|fixed);M23(Oxidation|15.994915|fixed)
        if mod_input != "":
            splitted_Modifications = mod_input.split(";")
            for mod in splitted_Modifications:

                (
                    position_or_aa_and_pos_unimod_name,
                    mod_mass,
                    fixed_or_opt,
                ) = mod.split("|")
                (
                    position_or_aa_and_pos,
                    unimod_name,
                ) = position_or_aa_and_pos_unimod_name.split("(")
                position_or_aa_and_pos = position_or_aa_and_pos.strip()
                unimod_name = unimod_name.strip()

                if position_or_aa_and_pos.upper() == "N-TERM":
                    position = 0
                else:
                    position = position_or_aa_and_pos[1:]

                translated_mods.append("{0}:{1}".format(unimod_name, position))

        modstring = ";".join(translated_mods)
        return modstring

    def get_column_names(self):
        """Get column names from param mapper.

        Returns:
            dict: dict mapping new to old column names
        """
        # create own uparma mapper
        headers = self.param_mapper.get_default_params(style=self.style)[
            "header_translations"
        ]["translated_value"]
        return headers
