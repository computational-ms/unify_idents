"""Engine parser."""
import xml.etree.ElementTree as ETree

import dask.dataframe as dd
import pandas as pd
import regex as re

from unify_idents.engine_parsers.base_parser import IdentBaseParser


def _get_single_spec_df(reference_dict, mapping_dict, spectrum):
    """Primary method for reading and storing information from a single spectrum.

    Args:
        reference_dict (dict): dict with reference columns to be filled in
        mapping_dict (dict): mapping of engine level column names to ursgal unified column names
        spectrum (xml Element): namespace of single spectrum with potentially multiple PSMs

    Returns:
        (pd.DataFrame): dataframe containing spectrum information

    """
    spec_records = []
    spec_level_dict = reference_dict.copy()
    spec_level_info = mapping_dict.keys() & spectrum.attrib.keys()
    spec_level_dict.update(
        {mapping_dict[k]: spectrum.attrib[k] for k in spec_level_info}
    )

    if "z" not in spectrum.attrib:
        return None
    spec_title = spectrum.findall('.//**[@label="Description"]')[0].text.split()[0]
    spec_level_dict["Spectrum Title"] = spec_title
    spec_level_dict["Spectrum ID"] = spec_title.split(".")[1]

    # Iterate children
    for psm in spectrum.findall(".//protein/*/domain"):
        psm_level_dict = spec_level_dict.copy()

        psm_level_dict["Calc m/z"] = psm.attrib["mh"]

        psm_level_info = mapping_dict.keys() & psm.attrib.keys()
        psm_level_dict.update({mapping_dict[k]: psm.attrib[k] for k in psm_level_info})

        # Record modifications
        mods = []
        for m in psm.findall(".//aa"):
            mass, abs_pos = m.attrib["modified"], m.attrib["at"]
            # abs pos is pos in protein, rel pos is pos in peptide
            rel_pos = int(abs_pos) - int(psm.attrib["start"])
            mods.append(f"{mass}:{rel_pos}")

        psm_level_dict["Modifications"] = mods

        spec_records.append(psm_level_dict)

    return pd.DataFrame(spec_records)


class XTandemAlanine_Parser(IdentBaseParser):
    """File parser for X!Tandem Alanine."""

    def __init__(self, *args, **kwargs):
        """Initialize parser.

        Reads in data file and provides mappings.
        """
        super().__init__(*args, **kwargs)
        self.style = "xtandem_style_1"
        tree = ETree.parse(self.input_file)
        self.root = tree.getroot()
        self.reference_dict["Search Engine"] = (
            "xtandem_"
            + re.search(
                r"(?<=Tandem )\w+",
                self.root.find(
                    './/*[@label="performance parameters"]/*[@label="process, version"]'
                ).text,
            )
            .group()
            .lower()
        )
        self.mapping_dict = {
            v: k
            for k, v in self.param_mapper.get_default_params(style=self.style)[
                "header_translations"
            ]["translated_value"].items()
        }
        self.reference_dict.update({k: None for k in self.mapping_dict.values()})

    @classmethod
    def check_parser_compatibility(cls, file):
        """Assert compatibility between file and parser.

        Args:
            file (str): path to input file

        Returns:
            bool: True if parser and file are compatible

        """
        is_xml = file.as_posix().endswith(".xml")
        with open(file.as_posix()) as f:
            try:
                head = "".join([next(f) for _ in range(10)])
            except StopIteration:
                head = ""
        contains_ref = "tandem-style.xsl" in head

        return is_xml and contains_ref

    def map_mod_names(self):
        """Map modification names in unify style.

        Args:
            df (pd.DataFrame): input dataframe
        """
        unique_mods = set().union(
            *self.df["Modifications"]
            .apply(set, meta=("Modifications", str))
            .compute()
            .values
        )
        unique_mod_masses = {m.split(":")[0] for m in unique_mods}
        potential_names = {
            m: [
                name
                for name in self.mod_mapper.mass_to_names(
                    round(float(m), 4), decimals=4
                )
                if name in self.mod_dict
            ]
            for m in unique_mod_masses
        }
        mod_translation = {}
        remap_dict = {}
        for m in unique_mods:
            mass, pos = m.split(":")
            potential_mods = potential_names[mass]
            if len(potential_mods) == 0:
                mod_translation[m] = None
            else:
                for name in potential_mods:
                    if ("Prot-N-term" in self.mod_dict[name]["position"]) and (
                        int(pos) == 0
                    ):
                        remap_dict[m] = f"{name}:0"
                    else:
                        remap_dict[m] = f"{name}:{int(pos)+1}"
        self.df["Modifications"] = self.df["Modifications"].apply(
            lambda x: ";".join(x), meta=("Modifications", str)
        )
        for old, new in remap_dict.items():
            self.df["Modifications"] = self.df["Modifications"].str.replace(
                old, new, regex=False
            )

    def unify(self):
        """
        Primary method to read and unify engine output.

        Returns:
            self.df (pd.DataFrame): unified dataframe
        """
        self.root = [element for element in self.root if "id" in element.keys()]
        data_structure = {
            ("df", i): (
                _get_single_spec_df,
                self.reference_dict,
                self.mapping_dict,
                spec,
            )
            for i, spec in enumerate(self.root)
        }
        df_type_mapping = [
            (k, self.dtype_mapping[k]) if k in self.dtype_mapping else (k, str)
            for k in self.reference_dict.keys()
        ]
        self.df = dd.DataFrame(
            data_structure, "df", df_type_mapping, (len(self.root) + 1) * [None]
        )
        self.df = self.df.repartition(partition_size="100MB")
        self.df["Calc m/z"] = (
            (self.df["Calc m/z"].astype(float) - self.PROTON)
            / self.df["Charge"].astype(int)
        ) + self.PROTON
        self.df = self.df.persist()
        self.map_mod_names()
        self.process_unify_style()

        return self.df
