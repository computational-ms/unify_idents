"""Parser handler."""
import multiprocessing as mp

import numpy as np
import pandas as pd
import regex as re
import uparma
from chemical_composition import ChemicalComposition
from loguru import logger
from peptide_mapper.mapper import UPeptideMapper
from unimod_mapper.unimod_mapper import UnimodMapper

from unify_idents.utils import merge_and_join_dicts
from itertools import repeat


def get_mass_and_composition(cc, seq, mods):
    """Compute theoretical mass and hill_notation of any single peptidoform.

    Returns None if sequence contains unknown amino acids.
    Args:
        cc (ChemicalComposition): chemical composition object instantiated with correct mods xml
        seq (str): peptide sequence
        mods (str): modifications of the peptide sequence, given as "UnimodName:Position"

    Returns:
        tuple: (computed mass, hill_notation_unimod string)

    """
    try:
        cc.use(sequence=seq, modifications=mods)
    except (KeyError, Exception):
        return None, None
    return cc.mass(), cc.hill_notation_unimod()


class BaseParser:
    """Base class of all parser types."""

    def __init__(self, input_file, params):
        """Initialize parser.

        Args:
            input_file (str): path to input file
            params (dict): ursgal param dict
        """
        self.input_file = input_file
        if params is None:
            params = {}
        self.params = params
        self.xml_file_list = self.params.get("xml_file_list", None)
        self.param_mapper = uparma.UParma()
        self.style = None

    @classmethod
    def check_parser_compatibility(cls, file):
        """Assert compatibility between file and parser.

        Args:
            file (str): path to input file

        Returns:
            bool: True if parser and file are compatible

        """
        return False


class IdentBaseParser(BaseParser):
    """Base class of all ident parsers."""

    def __init__(self, *args, **kwargs):
        """Initialize parser.

        Reads in data file and provides mappings.
        """
        super().__init__(*args, **kwargs)
        self.DELIMITER = self.params.get("delimiter", "<|>")
        self.PROTON = 1.00727646677
        self.df = None
        self.mod_mapper = UnimodMapper(xml_file_list=self.xml_file_list)
        self.params["mapped_mods"] = self.mod_mapper.map_mods(
            mod_list=self.params.get("modifications", [])
        )
        self.mod_dict = self._create_mod_dicts()
        self.rt_truncate_precision = 2
        self.reference_dict = {
            "exp_mz": None,
            "calc_mz": None,
            "spectrum_title": None,
            "search_engine": None,
            "spectrum_id": None,
            "modifications": None,
            "retention_time_seconds": None,
        }
        self.dtype_mapping = {
            "spectrum_title": "str",
            "raw_data_location": "str",
            "spectrum_id": "int32",
            "sequence": "str",
            "modifications": "str",
            "charge": "int32",
            "is_decoy": "bool",
            "rank": "int32",
            "protein_id": "str",
            "retention_time_seconds": "float32",
            "exp_mz": "float32",
            "calc_mz": "float32",
            "ucalc_mz": "float32",
            "ucalc_mass": "float32",
            "accuracy_ppm": "float32",
            "chemical_composition": "str",
            "sequence_start": "str",
            "sequence_stop": "str",
            "sequence_pre_aa": "str",
            "sequence_post_aa": "str",
            "enzn": "str",
            "enzc": "str",
            "missed_cleavages": "int32",
            "search_engine": "str",
        }
        self.col_order = pd.Series(self.dtype_mapping.keys())

    def _calc_mz(self, mass, charge):
        """Calulates mass-to-charge ratio.

        Args:
            mass (pd.Series): masses
            charge (pd.Series): charges

        Returns:
            (pd.Series): m/z
        """
        return (
            mass.astype(float) + (charge.astype(int) * self.PROTON)
        ) / charge.astype(int)

    def _create_mod_dicts(self):
        """
        Create dict containing meta information about static and variable mods.

        Returns:
            mod_dict (dict): mapped modifications and information
        """
        mod_dict = {}
        for mod_type in ["fix", "opt"]:
            for modification in self.params["mapped_mods"][mod_type]:
                aa = modification["aa"]
                pos = modification["position"]
                name = modification["name"]
                if name not in mod_dict.keys():
                    mod_dict[name] = {
                        "mass": modification["mass"],
                        "aa": set(),
                        "position": set(),
                    }
                mod_dict[name]["aa"].add(aa)

                mod_dict[name]["aa"].add(pos)
                mod_dict[name]["position"].add(pos)

        return mod_dict

    def clean_up_modifications(self):
        """Sanitizes modstrings generated by engine parsers.

        Modifications are sorted by position and leading, repeated or trailing delimiters are removed
        Operations are performed inplace on self.df
        """
        # Remove any trailing or leading delimiters or only-delimiter modstrings
        self.df.loc[:, "modifications"] = self.df.loc[:, "modifications"].str.replace(
            r"^;+(?=\w)", "", regex=True
        )
        self.df.loc[:, "modifications"] = self.df.loc[:, "modifications"].str.replace(
            r"(?<=\w);+$", "", regex=True
        )
        self.df.loc[:, "modifications"] = self.df.loc[:, "modifications"].str.replace(
            r"^;+$", "", regex=True
        )

        # Ensure same order of modifications

        self.df.loc[:, "modifications"] = (
            self.df["modifications"].fillna("").str.split(";").apply(self.sort_mods)
        )

    def sort_mods(self, data):
        """
        Sort mods based on position in peptide.

        Args:
            data (list): list of modifications with style "Mod:position"
        Returns:
            sorted_formatted_mods (str): String with sorted mods in style "Mod1:pos1;Modn:posn"
        """
        sort_pattern = r"(\w+)(?:\:)(\d+)"
        positions = [int(re.search(sort_pattern, d).group(2)) for d in data if d != ""]
        names = [re.search(sort_pattern, d).group(1) for d in data if d != ""]
        sorted_mods = sorted(zip(names, positions), key=lambda x: x[1])
        sorted_mods_str = [[str(i) for i in m] for m in sorted_mods]
        sorted_formatted_mods = ";".join([":".join(m) for m in sorted_mods_str])
        return sorted_formatted_mods

    def assert_only_iupac_and_missing_aas(self):
        """Assert that only IUPAC nomenclature one letter amino acids are used in sequence.

        Non-IUPAC designations are dropped.
        Operations are performed inplace.
        """
        self.df["sequence"] = self.df["sequence"].str.upper()
        # Added X for missing AAs
        iupac_aas = set("ACDEFGHIKLMNPQRSTVWY")
        iupac_conform_seqs = self.df["sequence"].apply(
            lambda seq: set(seq).issubset(iupac_aas)
        )
        if any(~iupac_conform_seqs):
            self.df = self.df.loc[iupac_conform_seqs, :]
            logger.warning(
                f"sequences are not IUPAC conform. {(~iupac_conform_seqs).sum()} PSMs were dropped."
            )

    def add_protein_ids(self):
        """Add all Protein IDs that matching the sequence.

        Operations are performed inplace on self.df
        """
        peptide_mapper = UPeptideMapper(self.params["database"])
        mapped_peptides = peptide_mapper.map_peptides(self.df["sequence"].tolist())

        peptide_mappings = [
            merge_and_join_dicts(mapped_peptides[seq], self.DELIMITER)
            for seq in self.df["sequence"]
        ]

        columns_translations = {
            "start": "sequence_start",
            "end": "sequence_stop",
            "post": "sequence_post_aa",
            "id": "protein_id",
            "pre": "sequence_pre_aa",
        }
        new_columns = pd.DataFrame(peptide_mappings)
        new_columns.rename(columns=columns_translations, inplace=True)

        self.df.loc[:, new_columns.columns] = new_columns.values
        new_columns = new_columns.dropna(axis=0, how="all")
        if len(new_columns) != len(self.df):
            logger.warning(
                f"{len(self.df)-len(new_columns)} PSMs were dropped because their respective sequences could not be mapped."
            )
        self.df = self.df.iloc[new_columns.index, :].reset_index(drop=True)

    def check_enzyme_specificity(self):
        """Check consistency of N/C-terminal cleavage sites.

        Calculates number of missed cleavage sites.
        Operations are performed inplace.
        """
        if self.params["enzyme"] == ".^":
            self.df.loc[:, ["enzn", "enzc"]] = True
            self.df.loc[:, "missed_cleavages"] = 0
            return None

        enzyme_pattern = self.params["enzyme"]
        integrity_strictness = self.params["terminal_cleavage_site_integrity"]

        pren_seq = (
            pd.concat(
                [
                    self.df["sequence_pre_aa"].str.split(rf"{self.DELIMITER}"),
                    self.df["sequence"].str[:1],
                ],
                axis=1,
            )
            .explode("sequence_pre_aa")
            .sum(axis=1)
        )
        self.df.loc[:, "enzn"] = (
            pren_seq.str.split(rf"{enzyme_pattern}").str[0].str.len() == 1
        ).groupby(pren_seq.index).agg(integrity_strictness) | (
            pren_seq.str[0] == "-"
        ).groupby(
            pren_seq.index
        ).agg(
            integrity_strictness
        )
        postc_seq = (
            pd.concat(
                [
                    self.df["sequence"].str[-1:],
                    self.df["sequence_post_aa"].str.split("<\\|>"),
                ],
                axis=1,
            )
            .explode("sequence_post_aa")
            .sum(axis=1)
        )
        self.df.loc[:, "enzc"] = (
            postc_seq.str.split(rf"{enzyme_pattern}").str[0].str.len() == 1
        ).groupby(postc_seq.index).agg(integrity_strictness) | (
            postc_seq.str[-1] == "-"
        ).groupby(
            postc_seq.index
        ).agg(
            integrity_strictness
        )

        internal_cuts = self.df["sequence"].str.split(rf"{enzyme_pattern}")
        self.df.loc[:, "missed_cleavages"] = (
            internal_cuts.apply(len)
            - internal_cuts.apply(lambda row: "" in row).astype(int)
            - 1
        )

    def calc_masses_offsets_and_composition(self):
        """Theoretical masses and mass-to-charge ratios are computed and added.

        Offsets are calculated between theoretical and experimental mass-to-charge ratio.
        Operations are performed inplace on self.df
        """
        with mp.Pool(self.params.get("cpus", mp.cpu_count() - 1)) as pool:
            cc_masses_and_comp = pool.starmap(
                get_mass_and_composition,
                zip(
                    repeat(
                        ChemicalComposition(
                            unimod_file_list=self.params.get("xml_file_list", None)
                        )
                    ),
                    self.df["sequence"].values,
                    self.df["modifications"].values,
                ),
                chunksize=1,
            )
        self.df.loc[:, ["ucalc_mass", "chemical_composition"]] = cc_masses_and_comp
        self.df.loc[:, "ucalc_mz"] = self._calc_mz(
            mass=self.df["ucalc_mass"], charge=self.df["charge"]
        )
        self.df.loc[:, "accuracy_ppm"] = (
            (self.df["exp_mz"].astype(float) - self.df["ucalc_mz"])
            / self.df["ucalc_mz"]
            * 1e6
        )

    def _read_meta_info_lookup_file(self):
        """Read meta info lookup file.

        Returns:
            rt_lookup (pd.DataFrame): loaded rt_pickle_file indexable by Spectrum ID
        """
        rt_lookup = pd.read_csv(self.params["rt_pickle_name"], compression="infer")
        rt_lookup["rt_unit"] = rt_lookup["rt_unit"].replace(
            {"second": 1, "minute": 60, "s": 1, "min": 60}
        )
        rt_lookup.set_index(
            [
                "spectrum_id",
                rt_lookup["rt"].apply(np.trunc, args=(self.rt_truncate_precision,)),
            ],
            inplace=True,
        )
        return rt_lookup

    def get_meta_info(self):
        """Extract meta information.

        Experimental mass-to-charge ratios, retention times, file names,
        and spectrum titles are added.
        Operations are performed inplace on self.df
        """
        rt_lookup = self._read_meta_info_lookup_file()
        spec_ids = self.df["spectrum_id"].astype(int)
        if self.style in ("comet_style_1", "omssa_style_1"):
            logger.warning(
                "This engine does not provide retention time information. Grouping only by Spectrum ID. This may cause problems when working with multi-file inputs."
            )
            rt_lookup = rt_lookup.loc[pd.IndexSlice[spec_ids.unique(), :]].droplevel(
                "rt"
            )
            spec_rt_idx = spec_ids
        else:
            spec_rt_idx = (
                pd.concat(
                    [
                        spec_ids,
                        self.df["retention_time_seconds"]
                        .astype(float)
                        .apply(np.trunc, args=(self.rt_truncate_precision,)),
                    ],
                    axis=1,
                )
                .apply(tuple, axis=1)
                .to_list()
            )
        try:
            self.df["retention_time_seconds"] = (
                rt_lookup.loc[spec_rt_idx, ["rt", "rt_unit"]]
                .astype(float)
                .product(axis=1)
                .to_list()
            )
        except KeyError:
            logger.warning("PSMs could not be uniquely mapped to meta information.")
            logger.info("Continuing with smallest delta retention times.")
            missing_truncated_indices = set(
                ind for ind in spec_rt_idx if ind not in rt_lookup.index
            )
            ind_mapping = {}
            for ind in missing_truncated_indices:
                meta_rt = rt_lookup.loc[ind[0], "rt"]
                smallest_delta_idx = abs(meta_rt - ind[1]).idxmin()
                if abs(meta_rt.loc[smallest_delta_idx] - ind[1]) <= 10 ** (
                    -self.rt_truncate_precision
                ):
                    ind_mapping[ind] = (ind[0], smallest_delta_idx)
                else:
                    logger.error(
                        f"No PSMs with (spectrum_id, retention_time) {ind} in meta information."
                    )
                    raise KeyError
            spec_rt_idx = [ind_mapping.get(ind, ind) for ind in spec_rt_idx]
            self.df["retention_time_seconds"] = (
                rt_lookup.loc[spec_rt_idx, ["rt", "rt_unit"]]
                .astype(float)
                .product(axis=1)
                .to_list()
            )

        self.df["exp_mz"] = rt_lookup.loc[spec_rt_idx, "precursor_mz"].to_list()
        self.df["raw_data_location"] = rt_lookup.loc[spec_rt_idx, "file"].to_list()
        self.df.loc[:, "spectrum_title"] = (
            self.df["raw_data_location"]
            + "."
            + self.df["spectrum_id"].astype(str)
            + "."
            + self.df["spectrum_id"].astype(str)
            + "."
            + self.df["charge"].astype(str)
        )

    def add_ranks(self):
        """Ranks are calculated based on the engine scoring column at Spectrum ID level.

        Operations are performed inplace on self.df
        """
        eng_name = self.df["search_engine"].unique()[0]
        score_col = self.params["validation_score_field"][eng_name]
        top_is_highest = self.params["bigger_scores_better"][eng_name]
        ranking_needs_to_be_ascending = False if top_is_highest is True else True
        self.df.loc[:, score_col] = self.df[score_col].astype(float)
        self.df.loc[:, "rank"] = (
            self.df.groupby("spectrum_id")[score_col]
            .rank(ascending=ranking_needs_to_be_ascending, method="min")
            .astype(int)
        )

    def add_decoy_identity(self):
        """Add boolean decoy state if designated decoy prefix is in Protein IDs.

        Operations are performed inplace on self.df
        """
        decoy_tag = self.params.get("decoy_tag", "decoy_")
        self.df.loc[:, "is_decoy"] = self.df["protein_id"].str.contains(decoy_tag)

    def sanitize(self):
        """Perform dataframe sanitation steps.

        - Missing raw data locations are filled with empty strings
        - Columns that were not filled in but should exist in the unified format are added and set to None
        - Columns in the dataframe which could not be properly mapped are removed (warning is raised)
        Operations are performed inplace on self.df
        """
        missing_data_locs = ~(self.df["raw_data_location"].str.len() > 0)
        self.df.loc[missing_data_locs, "raw_data_location"] = ""

        # Set missing columns to None and reorder columns in standardized manner
        new_cols = self.col_order[~self.col_order.isin(self.df.columns)].to_list()
        self.df.loc[:, new_cols] = None
        self.df = self.df.loc[
            :,
            self.col_order.tolist()
            + sorted(self.df.columns[~self.df.columns.isin(self.col_order)].tolist()),
        ]
        self.df = self.df.astype(self.dtype_mapping)

        # Ensure there are not any column that should not be
        if hasattr(self, "mapping_dict"):
            new_cols = set(self.mapping_dict.keys())
        else:
            new_cols = set()
        additional_cols = set(self.df.columns).difference(
            set(self.dtype_mapping.keys()) | new_cols
        )
        unmapped_add_cols = [c for c in additional_cols if ":" not in c]
        if len(unmapped_add_cols) > 0:
            logger.warning(
                f"Some engine level columns ({unmapped_add_cols}) were not properly mapped and removed."
            )
            self.df.drop(columns=unmapped_add_cols, inplace=True, errors="ignore")

        # Drop unwanted duplicated rows
        init_len = len(self.df)
        self.df.drop_duplicates(inplace=True)
        rows_dropped = init_len - len(self.df)
        if rows_dropped != 0:
            logger.warning(
                f"{rows_dropped} duplicated rows were dropped in output csv."
            )

    def process_unify_style(self):
        """Combine all additional operations that are needed to calculate new columns and sanitize the dataframe.

        Operations are performed inplace on self.df
        """
        self.df.drop_duplicates(inplace=True, ignore_index=True)
        self.clean_up_modifications()
        self.assert_only_iupac_and_missing_aas()
        self.add_protein_ids()
        self.get_meta_info()
        self.calc_masses_offsets_and_composition()
        self.check_enzyme_specificity()
        self.add_ranks()
        self.add_decoy_identity()
        self.sanitize()


class QuantBaseParser(BaseParser):
    """Base class of all quant parsers."""

    def __init__(self, *args, **kwargs):
        """Initialize parser.

        Reads in data file and provides mappings.
        """
        super().__init__(*args, **kwargs)
        self.cc = ChemicalComposition()
        self.required_headers = {
            "file_name",
            "spectrum_id",
            "trivial_name",
            "chemical_composition",
            "precursor_spectrum_id",
            "retention_time",
            "charge",
            "quant_run_id",
            "quant_value",
            "quant_score",
            "quant_group",
            "processing_level",
            "delta_mz",
            "label",
            "condition",
            "ident_reference",
            "fwhm",
            "s2i",
            "p2t",
            "coalescence",
        }

    def process_unify_style(self):
        """Apply sanitizing methods."""
        pass
