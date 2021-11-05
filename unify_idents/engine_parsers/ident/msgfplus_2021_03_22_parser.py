import xml.etree.ElementTree as ETree
from pathlib import Path
import regex as re
import multiprocessing as mp
import pandas as pd
from tqdm import tqdm
from itertools import repeat

from unify_idents.engine_parsers.base_parser import __IdentBaseParser


def _get_single_spec_df(reference_dict, spectrum):
    spec_records = []
    spec_level_dict = reference_dict.copy()
    spec_level_info = {
        i.attrib["name"]: i.attrib["value"]
        for i in spectrum.findall(".//{*}cvParam")
        if i.attrib["name"] in ["scan number(s)", "spectrum title"]
    }
    spec_level_dict["Spectrum ID"] = spec_level_info["scan number(s)"]
    spec_level_dict["Spectrum Title"] = spec_level_info["spectrum title"]

    # Iterate children
    for psm in spectrum.findall(".//{*}SpectrumIdentificationItem"):
        psm_level_dict = spec_level_dict.copy()

        # peptide_info = peptide_lookup[psm.attrib["peptide_ref"]]

        psm_level_dict.update(
            {
                "Exp m/z": psm.attrib["experimentalMassToCharge"],
                "Calc m/z": psm.attrib["calculatedMassToCharge"],
                "Charge": psm.attrib["chargeState"],
                "Sequence": psm.attrib["peptide_ref"],
                # "Sequence": peptide_info["Sequence"],
                # "Modifications": ";".join(peptide_info["Modifications"])
            }
        )
        psm_level_dict.update(
            {c.attrib["name"]: c.attrib["value"] for c in psm.findall(".//{*}cvParam")}
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
        cols_to_remove = [
            # TODO: this shouldnt exist in uparma or is mapping just wrong?
            "proteinacc_start_stop_pre_post_;",
            # TODO: should this be in here?
            "Raw data location",
        ]
        self.mapping_dict = {
            v: k
            for k, v in self.param_mapper.get_default_params(style=self.style)[
                "header_translations"
            ]["translated_value"].items()
            if k not in cols_to_remove
        }
        self.reference_dict.update({k: None for k in self.mapping_dict.values()})

    @classmethod
    def check_parser_compatibility(cls, file):
        is_mzid = file.as_posix().endswith(".mzid")

        with open(file.as_posix()) as f:
            head = "".join([next(f) for x in range(20)])
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
                    zip(repeat(self.reference_dict), spec_idents),
                    total=len(spec_idents),
                ),
                chunksize=1,
            )
        unified_df = pd.concat(chunk_dfs, axis=0, ignore_index=True)
        seq_mods = pd.DataFrame(unified_df["Sequence"].map(peptide_lookup).to_list())
        unified_df.loc[:, seq_mods.columns] = seq_mods

        return unified_df
