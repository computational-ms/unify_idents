"""Engine parser."""
import multiprocessing as mp
import sys
import xml.etree.ElementTree as ETree
from itertools import repeat

import pandas as pd
from loguru import logger
from tqdm import tqdm

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
        self.reference_dict.update(
            {
                "Raw data location": self.root.attrib["label"]
                .split("models from ")[1]
                .replace("'", ""),
                "Search Engine": "xtandem_alanine",
            }
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
            head = "".join([next(f) for _ in range(10)])
        contains_ref = "tandem-style.xsl" in head

        return is_xml and contains_ref

    def map_mod_names(self, df):
        """Map modification names in unify style.

        Args:
            df (pd.DataFrame): input dataframe

        Returns:
            df (pd.DataFrame): dataframe with processed modification column

        """
        unique_mods = set().union(*df["Modifications"].apply(set).values)
        unique_mod_masses = {m.split(":")[0] for m in unique_mods}
        potential_names = {
            m: [
                name
                for name in self.mod_mapper.mass_to_names(round(float(m), 4), decimals=4)
                if name in self.mod_dict
            ]
            for m in unique_mod_masses
        }
        mod_translation = {}
        new_mods = pd.Series("", index=df.index)
        for m in unique_mods:
            mass, pos = m.split(":")
            potential_mods = potential_names[mass]
            if len(potential_mods) == 0:
                mod_translation[m] = None
            else:
                for name in potential_mods:
                    # TODO: Is position 'any' respected here
                    in_seq = df["Sequence"].str[int(pos)].isin(
                        self.mod_dict[name]["aa"]
                    ) & df["Modifications"].str.join("|").str.contains(m)
                    if in_seq.sum() != 0:
                        new_mods.loc[in_seq] += f"{name}:{int(pos)+1};"
                    n_term = (~in_seq) & (
                        ("Prot-N-term" in self.mod_dict[name]["position"])
                        & df["Modifications"].str.join("|").str.contains(m)
                    )
                    if n_term.sum() != 0:
                        new_mods.loc[n_term] += f"{name}:0;"
        df["Modifications"] = new_mods.str.rstrip(";")

        return df

    def unify(self):
        """
        Primary method to read and unify engine output.

        Returns:
            self.df (pd.DataFrame): unified dataframe
        """
        logger.remove()
        logger.add(lambda msg: tqdm.write(msg, end=""))
        pbar_iterator = tqdm(
            zip(
                repeat(self.reference_dict),
                repeat(self.mapping_dict),
                self.root,
            ),
            total=len(self.root),
        )
        with mp.Pool(self.params.get("cpus", mp.cpu_count() - 1)) as pool:
            chunk_dfs = pool.starmap(_get_single_spec_df, pbar_iterator)
        logger.remove()
        logger.add(sys.stdout)
        chunk_dfs = [df for df in chunk_dfs if not df is None]
        self.df = pd.concat(chunk_dfs, axis=0, ignore_index=True)
        self.df["Calc m/z"] = (
            (self.df["Calc m/z"].astype(float) - self.PROTON)
            / self.df["Charge"].astype(int)
        ) + self.PROTON
        self.df = self.map_mod_names(self.df)
        self.process_unify_style()

        return self.df
