import multiprocessing as mp
import xml.etree.ElementTree as ETree
from itertools import repeat

import pandas as pd
import regex as re
from tqdm import tqdm

from unify_idents.engine_parsers.base_parser import __IdentBaseParser


def _get_single_spec_df(reference_dict, mapping_dict, spectrum):
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


class MSGFPlus_2021_03_22(__IdentBaseParser):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.style = "msgfplus_style_1"

        tree = ETree.parse(self.input_file)
        self.root = tree.getroot()
        self.reference_dict.update(
            {
                "Raw data location": self.root.find(".//{*}SpectraData").attrib[
                    "location"
                ],
                "Search Engine": "msgfplus_"
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
        is_mzid = file.as_posix().endswith(".mzid")

        with open(file.as_posix()) as f:
            try:
                head = "".join([next(f) for x in range(20)])
            except StopIteration:
                head = ""
        contains_engine = "MS-GF+" in head
        contains_correct_version = (
            "Release (v2021.03.22)" in head or "Release (v2019.07.03)" in head
        )
        return is_mzid and contains_engine and contains_correct_version

    def _get_peptide_lookup(self):
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
        peptide_lookup = self._get_peptide_lookup()
        spec_idents = self.root.findall(".//{*}SpectrumIdentificationResult")
        with mp.Pool(self.params.get("cpus", mp.cpu_count() - 1)) as pool:
            chunk_dfs = pool.starmap(
                _get_single_spec_df,
                tqdm(
                    zip(
                        repeat(self.reference_dict),
                        repeat(self.mapping_dict),
                        spec_idents,
                    ),
                    total=len(spec_idents),
                ),
                chunksize=1,
            )
        self.df = pd.concat(chunk_dfs, axis=0, ignore_index=True)
        seq_mods = pd.DataFrame(self.df["Sequence"].map(peptide_lookup).to_list())
        self.df.loc[:, seq_mods.columns] = seq_mods
        self.process_unify_style()

        return self.df
