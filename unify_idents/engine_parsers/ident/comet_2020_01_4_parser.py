"""Engine parser."""
import multiprocessing as mp
import sys
import xml.etree.ElementTree as ETree
from itertools import repeat

import numpy as np
import pandas as pd
import regex as re
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
    spec_level_dict["Spectrum ID"] = spectrum.attrib["spectrumID"].split("scan=")[-1]

    # Iterate children
    for psm in spectrum.findall(".//{*}SpectrumIdentificationItem"):
        psm_level_dict = spec_level_dict.copy()
        psm_level_dict.update(
            {mapping_dict[k]: psm.attrib[k] for k in mapping_dict if k in psm.attrib}
        )
        cv_param_info = {
            c.attrib["name"]: c.attrib["value"] for c in psm.findall(".//{*}cvParam")
        }
        psm_level_dict.update(
            {
                mapping_dict[k]: cv_param_info[k]
                for k in mapping_dict
                if k in cv_param_info
            }
        )

        spec_records.append(psm_level_dict)
    return pd.DataFrame(spec_records)


class Comet_2020_01_4_Parser(IdentBaseParser):
    """File parser for Comet."""

    def __init__(self, *args, **kwargs):
        """Initialize parser.

        Reads in data file and provides mappings.
        """
        super().__init__(*args, **kwargs)
        self.style = "comet_style_1"

        tree = ETree.parse(self.input_file)
        self.root = tree.getroot()
        self.reference_dict.update(
            {
                "Raw data location": self.root.find(
                    ".//{*}DataCollection/{*}Inputs/{*}SpectraData"
                ).attrib["location"],
                "Search Engine": "comet_"
                + "_".join(
                    re.findall(
                        r"([/d]*\d+)",
                        self.root.find(".//{*}AnalysisSoftware").attrib["version"],
                    )
                ),
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
        is_mzid = file.as_posix().endswith(".mzid")

        with open(file.as_posix()) as f:
            try:
                head = "".join([next(f) for _ in range(10)])
            except StopIteration:
                head = ""
        contains_engine = "Comet" in head
        return is_mzid and contains_engine

    def _map_mods_and_sequences(self):
        """Replace internal tags to retrieve sequences and formatted modification strings.

        Operations are performed inplace.
        """
        # Register fixed mods
        modifications = self.root.findall(
            ".//{*}AnalysisProtocolCollection/{*}SpectrumIdentificationProtocol/{*}ModificationParams/{*}SearchModification"
        )
        fixed_mods = {
            sm.attrib["residues"]: sm.find(".//{*}cvParam[@cvRef='UNIMOD']").attrib[
                "name"
            ]
            for sm in modifications
            if sm.attrib["fixedMod"] == "true"
        }
        if len(fixed_mods) > 0:
            fixed_mod_strings = []
            for fm_res, fm_name in fixed_mods.items():
                fixed_mod_strings.append(
                    self.df["Sequence"]
                    .str.split(fm_res)
                    .apply(
                        lambda l: ";".join(
                            [
                                fm_name + ":" + ind
                                for ind in (
                                    np.cumsum(list(map(len, l[:-1]))) + range(1, len(l))
                                ).astype(str)
                            ]
                        )
                    )
                )

            fixed_mod_strings = (
                pd.concat(fixed_mod_strings, axis=1)
                .agg(";".join, axis=1)
                .str.rstrip(";")
            )

        modification_mass_map = {
            sm.attrib["massDelta"]: sm.find(".//{*}cvParam[@cvRef='UNIMOD']").attrib[
                "name"
            ]
            for sm in modifications
        }
        lookup = {}
        for pep in self.root.findall(".//{*}Peptide"):
            id = pep.attrib.get("id", "")
            lookup[id] = {"Modifications": []}
            lookup[id]["Sequence"] = pep.find(".//{*}PeptideSequence").text
            for child in pep.findall(".//{*}Modification"):
                lookup[id]["Modifications"].append(
                    f"{modification_mass_map[child.attrib['monoisotopicMassDelta']]}:{child.attrib['location']}"
                )
            lookup[id]["Modifications"] = ";".join(lookup[id]["Modifications"])

        seq_mods = pd.DataFrame(self.df["Sequence"].map(lookup).to_list())
        self.df.loc[:, "Modifications"] = (
            seq_mods["Modifications"].str.cat(fixed_mod_strings, sep=";").str.strip(";")
        )
        self.df.loc[:, "Sequence"] = seq_mods["Sequence"]

    def unify(self):
        """
        Primary method to read and unify engine output.

        Returns:
            self.df (pd.DataFrame): unified dataframe
        """
        spec_idents = self.root.findall(
            ".//{*}SpectrumIdentificationList/{*}SpectrumIdentificationResult"
        )
        logger.remove()
        logger.add(lambda msg: tqdm.write(msg, end=""))
        pbar_iterator = tqdm(
            zip(
                repeat(self.reference_dict),
                repeat(self.mapping_dict),
                spec_idents,
            ),
            total=len(spec_idents),
        )
        with mp.Pool(self.params.get("cpus", mp.cpu_count() - 1)) as pool:
            chunk_dfs = pool.starmap(
                _get_single_spec_df,
                pbar_iterator,
            )
        logger.remove()
        logger.add(sys.stdout)
        self.df = pd.concat(chunk_dfs, axis=0, ignore_index=True)
        self._map_mods_and_sequences()
        self.df.loc[:, "Spectrum Title"] = (
            self.df["Raw data location"].str.extract(f"(?<=/)([\w_]+)(?=\.)")[0]
            + "."
            + self.df["Spectrum ID"]
            + "."
            + self.df["Spectrum ID"]
            + "."
            + self.df["Charge"]
        )
        self.process_unify_style()

        return self.df
