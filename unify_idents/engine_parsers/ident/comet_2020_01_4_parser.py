import multiprocessing as mp
import xml.etree.ElementTree as ETree
from itertools import repeat

import pandas as pd
import regex as re
from tqdm import tqdm

from unify_idents.engine_parsers.base_parser import __IdentBaseParser


def _get_single_spec_df(reference_dict, spectrum):
    spec_records = []
    spec_level_dict = reference_dict.copy()
    spec_level_dict["Spectrum ID"] = spectrum.attrib["spectrumID"].split("scan=")[-1]

    # Iterate children
    for psm in spectrum.findall(".//{*}SpectrumIdentificationItem"):
        psm_level_dict = spec_level_dict.copy()
        psm_level_dict.update(
            {
                "Exp m/z": psm.attrib["experimentalMassToCharge"],
                "Calc m/z": psm.attrib["calculatedMassToCharge"],
                "Charge": psm.attrib["chargeState"],
                "Sequence": psm.attrib["peptide_ref"],
            }
        )
        psm_level_dict.update(
            {c.attrib["name"]: c.attrib["value"] for c in psm.findall(".//{*}cvParam")}
        )

        spec_records.append(psm_level_dict)
    return pd.DataFrame(spec_records)


class Comet_2020_01_4_Parser(__IdentBaseParser):
    def __init__(self, *args, **kwargs):
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
        # cols_to_remove = [
        #     # TODO: this shouldnt exist in uparma or is mapping just wrong?
        #     "proteinacc_start_stop_pre_post_;",
        #     # TODO: should this be in here?
        #     "Raw data location",
        # ]
        # self.mapping_dict = {
        #     v: k
        #     for k, v in self.param_mapper.get_default_params(style=self.style)[
        #         "header_translations"
        #     ]["translated_value"].items()
        #     if k not in cols_to_remove
        # }
        # self.reference_dict.update({k: None for k in self.mapping_dict.values()})

    @classmethod
    def check_parser_compatibility(cls, file):
        is_mzid = file.as_posix().endswith(".mzid")

        with open(file.as_posix()) as f:
            try:
                head = "".join([next(f) for x in range(10)])
            except StopIteration:
                head = ""
        contains_engine = "Comet" in head
        contains_correct_version = "2020.01 rev. 4" in head
        return is_mzid and contains_engine and contains_correct_version

    def _get_peptide_lookup(self):
        modifications = self.root.findall(
            ".//{*}AnalysisProtocolCollection/{*}SpectrumIdentificationProtocol/{*}ModificationParams/{*}SearchModification"
        )
        modification_mass_map = {
            m.attrib["massDelta"]: m.find(".//{*}cvParam[@cvRef='UNIMOD']").attrib[
                "name"
            ]
            for m in modifications
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
        return lookup

    def unify(self):
        peptide_lookup = self._get_peptide_lookup()
        spec_idents = self.root.findall(
            ".//{*}SpectrumIdentificationList/{*}SpectrumIdentificationResult"
        )
        with mp.Pool(self.params.get("cpus", mp.cpu_count() - 1)) as pool:
            chunk_dfs = pool.starmap(
                _get_single_spec_df,
                tqdm(
                    zip(repeat(self.reference_dict), spec_idents),
                    total=len(spec_idents),
                ),
            )
        self.df = pd.concat(chunk_dfs, axis=0, ignore_index=True)
        seq_mods = pd.DataFrame(self.df["Sequence"].map(peptide_lookup).to_list())
        self.df.loc[:, seq_mods.columns] = seq_mods
        # TODO: what to do with the spectrum titles?
        self.process_unify_style()

        return self.df
