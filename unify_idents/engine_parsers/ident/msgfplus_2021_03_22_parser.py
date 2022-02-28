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
    spec_level_info = {
        c.attrib["name"]: c.attrib["value"] for c in spectrum.findall(".//{*}cvParam")
    }
    spec_level_dict.update(
        {
            mapping_dict[k]: spec_level_info[k]
            for k in mapping_dict
            if k in spec_level_info
        }
    )

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
        user_param_info = {
            c.attrib["name"]: c.attrib["value"] for c in psm.findall(".//{*}userParam")
        }
        psm_level_dict.update(
            {
                mapping_dict[k]: user_param_info[k]
                for k in mapping_dict
                if k in user_param_info
            }
        )

        spec_records.append(psm_level_dict)
    return pd.DataFrame(spec_records)


class MSGFPlus_2021_03_22_Parser(IdentBaseParser):
    """File parser for MSGF+."""

    def __init__(self, *args, **kwargs):
        """Initialize parser.

        Reads in data file and provides mappings.
        """
        super().__init__(*args, **kwargs)
        self.style = "msgfplus_style_1"

        tree = ETree.parse(self.input_file)
        self.root = tree.getroot()
        self.reference_dict["Search Engine"] = "msgfplus_" + "_".join(
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
                head = "".join([next(f) for _ in range(20)])
            except StopIteration:
                head = ""
        contains_engine = "MS-GF+" in head
        return is_mzid and contains_engine

    def _get_peptide_lookup(self):
        """Replace internal tags to retrieve sequences and formatted modification strings.

        Operations are performed inplace.
        """
        lookup = {}
        for pep in self.root.findall(".//{*}Peptide"):
            id = pep.attrib.get("id", "")
            lookup[id] = {"Modifications": []}
            for child in pep.findall(".//{*}PeptideSequence"):
                lookup[id]["Sequence"] = child.text
            for child in pep.findall(".//{*}Modification"):
                lookup[id]["Modifications"].append(
                    f"{child.find('.//{*}cvParam').attrib['name']}:{child.attrib['location']}"
                )
            lookup[id]["Modifications"] = ";".join(lookup[id]["Modifications"])
        return lookup

    def unify(self):
        """
        Primary method to read and unify engine output.

        Returns:
            self.df (pd.DataFrame): unified dataframe
        """
        peptide_lookup = self._get_peptide_lookup()
        spec_idents = self.root.findall(".//{*}SpectrumIdentificationResult")
        data_structure = {
            ("df", i): (
                _get_single_spec_df,
                self.reference_dict,
                self.mapping_dict,
                spec,
            )
            for i, spec in enumerate(spec_idents)
        }
        df_type_mapping = [
            (k, self.dtype_mapping[k]) if k in self.dtype_mapping else (k, str)
            for k in self.reference_dict.keys()
        ]
        self.df = dd.DataFrame(
            data_structure, "df", df_type_mapping, (len(spec_idents) + 1) * [None]
        )
        self.df = self.df.repartition(partition_size="100MB")
        seq_mods = (
            self.df["Sequence"]
            .map(peptide_lookup)
            .apply(pd.Series, meta=[("Modifications", str), ("Sequence", str)])
        )
        self.df[seq_mods.columns.to_list()] = seq_mods
        self.process_unify_style()

        return self.df
