import xml.etree.ElementTree as ETree

import numpy as np
import pandas as pd

from unify_idents.engine_parsers.base_parser import __IdentBaseParser


class OmssaParser(__IdentBaseParser):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.style = "omssa_style_1"

        self.df = pd.read_csv(self.input_file)

        cols_to_remove = [
            # TODO: this shouldnt exist in uparma or is mapping just wrong?
            "proteinacc_start_stop_pre_post_;",
            "Defline",
            # end todo
            "Start",
            "Stop",
            "NIST score",
            "gi",
            "Accession",
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
        self.df.drop(columns=cols_to_remove, inplace=True, errors="ignore")
        self.reference_dict.update({k: None for k in self.mapping_dict.values()})

    @classmethod
    def check_parser_compatibility(cls, file):
        is_csv = file.as_posix().endswith(".csv")
        with open(file.as_posix()) as f:
            head = "".join([next(f) for x in range(1)])
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

    def translate_mods(self):
        # TODO: Need sanitation already???
        self.df["Sequence"] = self.df["Sequence"].str.upper()
        fix_mods = []
        # Map fixed mods
        fixed_mod_types = [
            d for d in self.params["modifications"] if d["type"] == "fix"
        ]
        for fm in fixed_mod_types:
            fix_mods.append(
                self.df["Sequence"]
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

        unique_mods = (
            set()
            .union(
                *self.df["Modifications"]
                .fillna("")
                .str.split(" ,")
                .apply(lambda l: [x.split(":")[0] for x in l])
                .map(set)
            )
            .difference({""})
        )

        mod_translations = {
            k: {
                "unimod_id": None,
                "omssa_unimod_id": None,
                "unimod_name": None,
                "omssa_name": None,
                "aa_targets": [],
            }
            for k in unique_mods
        }
        for file in ["mods.xml", "usermods.xml"]:
            tree = ETree.parse(self.params["omssa_mod_dir"] / file)
            root = tree.getroot()
            for mod in unique_mods:
                # Get all parents of children where full text of a tag matches the modification str
                for element in root.findall(f".//*[.='{mod}']/.."):
                    for targ in element.findall(".//*"):
                        # Find names, match unimod id, record target aa
                        if "MSModSpec_unimod" in targ.tag:
                            mod_translations[mod]["omssa_unimod_id"] = targ.text
                        elif "MSModSpec_psi-ms" in targ.tag:
                            mod_translations[mod]["unimod_name"] = targ.text
                        elif "MSModSpec_residues" in targ.tag:
                            mod_translations[mod]["aa_targets"].extend(
                                [
                                    aa.text
                                    for aa in targ.findall(".//{*}MSModSpec_residues_E")
                                ]
                            )
                        elif "MSModSpec_name" in targ.tag:
                            mod_translations[mod]["omssa_name"] = targ.text
                            terminal_mod = ""
                            if "protein" in targ.tag:
                                terminal_mod += "Prot-"
                            if "n-term" in targ.tag:
                                terminal_mod += "N-term"
                            elif "c-term" in targ.tag:
                                terminal_mod += "C-term"
                            if len(terminal_mod) > 0:
                                mod_translations[mod]["aa_targets"].extend(
                                    [terminal_mod]
                                )

        # Apply fixes
        omssa_mod_corrections = {
            # this dict holds corrections of wrong OMSSA to unimod assignments
            "TMT 6-plex on K": {
                "unimod_id": "737",
                "omssa_unimod_id": "738",  # this is TMT duplex in unimod
                "unimod_name": "TMT6plex",
            },
            "TMT 6-plex on n-term peptide": {
                "unimod_id": "737",
                "omssa_unimod_id": "738",  # this is TMT duplex in unimod
                "unimod_name": "TMT6plex",
                "aa_targets": ["N-term"],  # override 'X' in OMSSA mods xml
            },
            "TMT duplex on K": {
                "unimod_id": "738",
                "omssa_unimod_id": "738",  # this is TMT duplex in unimod
                "unimod_name": "TMT2plex",
            },
            "TMT duplex on n-term peptide": {
                "unimod_id": "738",
                "omssa_unimod_id": "738",  # this is TMT duplex in unimod
                "unimod_name": "TMT2plex",
                "aa_targets": ["N-term"],  # override 'X' in OMSSA mods xml
            },
        }
        for k in mod_translations.keys():
            if k in omssa_mod_corrections.keys():
                mod_translations[k].update(omssa_mod_corrections[k])
            if mod_translations[k]["unimod_id"] is not None:
                mod_translations[k]["unimod_name"] = self.mod_mapper.id2name_list(
                    mod_translations[k]["unimod_id"]
                )[0]

        def _replace_mod_strings(row):
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

        opt_mods = (
            self.df["Modifications"].fillna("").apply(lambda r: _replace_mod_strings(r))
        )
        if len(fix_mods) > 0:
            fix_mods = pd.concat(fix_mods, axis=1)
            comb = pd.concat([opt_mods, fix_mods], axis=1)
            comb = (
                comb["Modifications"]
                .str.cat(comb["Sequence"], sep=";")
                .str.rstrip(";")
                .to_list()
            )
        else:
            comb = opt_mods.to_list()
        return comb

    def unify(self):
        self.df["Calc m/z"] = self._calc_mz(
            mass=self.df["Calc m/z"], charge=self.df["Charge"]
        )
        self.df["Spectrum ID"] = (
            self.df["Spectrum Title"].str.split(".").str[1].astype(int)
        )
        self.df["Raw data location"] = self.params.get(
            "Raw data location", self.df["Spectrum Title"].str.split(".").str[0]
        )
        self.df["Search Engine"] = "omssa_2_1_9"
        self.df["Modifications"] = self.translate_mods()
        self.process_unify_style()

        return self.df
