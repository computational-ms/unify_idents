"""Engine parser."""
import multiprocessing as mp
from io import BytesIO

import numpy as np
import pandas as pd
import regex as re
import sys
from loguru import logger
from lxml import etree
from tqdm import tqdm

from unify_idents.engine_parsers.ident.ident_base_parser import IdentBaseParser


def _mp_specs_init(func, reference_dict, mapping_dict):
    func.reference_dict = reference_dict
    func.mapping_dict = mapping_dict


def _get_single_spec_df(spectrum):
    """Primary method for reading and storing information from a single spectrum.

    Attributes:
        reference_dict (dict): dict with reference columns to be filled in
        mapping_dict (dict): mapping of engine level column names to ursgal unified column names

    Args:
        spectrum (xml Element): namespace of single spectrum with potentially multiple PSMs

    Returns:
        (pd.DataFrame): dataframe containing spectrum information

    """
    spectrum = etree.parse(BytesIO(spectrum)).getroot()
    spec_records = []
    spec_level_dict = _get_single_spec_df.reference_dict.copy()
    spec_level_dict["spectrum_id"] = spectrum.attrib["spectrumID"].split("scan=")[-1]

    # Iterate children
    for psm in spectrum.findall(".//{*}SpectrumIdentificationItem"):
        psm_level_dict = spec_level_dict.copy()
        psm_level_dict.update(
            {
                _get_single_spec_df.mapping_dict[k]: psm.attrib[k]
                for k in _get_single_spec_df.mapping_dict
                if k in psm.attrib
            }
        )
        cv_param_info = {
            c.attrib["name"]: c.attrib["value"] for c in psm.findall(".//{*}cvParam")
        }
        psm_level_dict.update(
            {
                _get_single_spec_df.mapping_dict[k]: cv_param_info[k]
                for k in _get_single_spec_df.mapping_dict
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

        tree = etree.parse(self.input_file)
        self.root = tree.getroot()
        self.reference_dict["search_engine"] = "comet_" + "_".join(
            re.findall(
                r"([/d]*\d+)",
                self.root.find(".//{*}AnalysisSoftware").attrib["version"],
            )
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
                    self.df["sequence"]
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
            lookup[id] = {"modifications": []}
            lookup[id]["sequence"] = pep.find(".//{*}PeptideSequence").text
            for child in pep.findall(".//{*}Modification"):
                lookup[id]["modifications"].append(
                    f"{modification_mass_map[child.attrib['monoisotopicMassDelta']]}:{child.attrib['location']}"
                )
            lookup[id]["modifications"] = ";".join(lookup[id]["modifications"])

        # TODO: check mod left strip
        seq_mods = pd.DataFrame(self.df["sequence"].map(lookup).to_list())
        self.df.loc[:, "modifications"] = (
            seq_mods["modifications"].str.cat(fixed_mod_strings, sep=";").str.strip(";")
        )
        self.df.loc[:, "sequence"] = seq_mods["sequence"]

    def unify(self):
        """
        Primary method to read and unify engine output.

        Returns:
            self.df (pd.DataFrame): unified dataframe
        """
        spec_idents = [
            etree.tostring(e)
            for e in self.root.findall(".//{*}SpectrumIdentificationResult")
        ]
        logger.remove()
        logger.add(lambda msg: tqdm.write(msg, end=""))
        with mp.Pool(
            self.params.get("cpus", mp.cpu_count() - 1),
            initializer=_mp_specs_init,
            initargs=(_get_single_spec_df, self.reference_dict, self.mapping_dict),
        ) as pool:
            chunk_dfs = pool.map(
                _get_single_spec_df,
                tqdm(spec_idents),
            )
        logger.remove()
        logger.add(sys.stdout)
        self.df = pd.concat(chunk_dfs, axis=0, ignore_index=True)
        self._map_mods_and_sequences()
        self.process_unify_style()

        return self.df
