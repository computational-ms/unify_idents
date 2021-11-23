"""Parser handler."""
import multiprocessing as mp

import pandas as pd
import uparma
from chemical_composition import ChemicalComposition
from loguru import logger
from peptide_mapper.mapper import UPeptideMapper
from unimod_mapper.unimod_mapper import UnimodMapper

from unify_idents.utils import merge_and_join_dicts

cc = ChemicalComposition()


def get_mass_and_composition(seq, mods):
    """Compute theoretical mass and hill_notation of any single peptidoform.

    Returns None if sequence contains unknown amino acids.
    Args:
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
        self.param_mapper = uparma.UParma()
        self.translated_params = self.param_mapper.get_default_params(
            "unify_csv_style_1"
        )
        self.translated_params.update(
            self.param_mapper.convert(self.params, "unify_csv_style_1")
        )

    @classmethod
    def check_parser_compatibility(cls, file):
        """Assert compatibility between file and parser.

        Args:
            file (str): path to input file

        Returns:
            bool: True if parser and file are compatible

        """
        return False

    def _read_rt_lookup_file(self):
        """Read retention time lookup file.

        Returns:
            rt_lookup (pd.DataFrame): loaded rt_pickle_file indexable by Spectrum ID
        """
        rt_lookup = pd.read_csv(self.params["rt_pickle_name"], compression="bz2")
        rt_lookup.set_index("Spectrum ID", inplace=True)
        rt_lookup["Unit"] = rt_lookup["Unit"].replace({"second": 1, "minute": 60})
        return rt_lookup


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
        self.mod_mapper = UnimodMapper()
        self.params["mapped_mods"] = self.mod_mapper.map_mods(
            mod_list=self.params.get("modifications", [])
        )
        self.mod_dict = self._create_mod_dicts()
        self.reference_dict = {
            "Exp m/z": None,
            "Calc m/z": None,
            "Spectrum Title": None,
            "Raw data location": None,
            "Search Engine": None,
            "Spectrum ID": None,
            "Modifications": None,
            "Retention Time (s)": None,
        }
        self.dtype_mapping = {
            "Spectrum Title": "str",
            "Raw data location": "str",
            "Spectrum ID": "int32",
            "Sequence": "str",
            "Modifications": "str",
            "Charge": "int32",
            "Is decoy": "bool",
            "Rank": "int32",
            "Protein ID": "str",
            "Retention Time (s)": "float32",
            "Exp m/z": "float32",
            "Calc m/z": "float32",
            "uCalc m/z": "float32",
            "uCalc Mass": "float32",
            "Accuracy (ppm)": "float32",
            "Mass Difference": "float32",
            "Chemical Composition": "str",
            "Sequence Start": "str",
            "Sequence Stop": "str",
            "Sequence Pre AA": "str",
            "Sequence Post AA": "str",
            "enzN": "str",
            "enzC": "str",
            "Missed Cleavages": "int32",
            "Search Engine": "str",
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
        return (mass.astype(float) + (charge.astype(int) * self.PROTON)) / charge.astype(
            int
        )

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

    def assert_only_iupac_and_missing_aas(self):
        """Assert that only IUPAC nomenclature one letter amino acids are used in sequence.

        Non-IUPAC designations are dropped.
        Operations are performed inplace.
        """
        self.df["Sequence"] = self.df["Sequence"].str.upper()
        # Added X for missing AAs
        iupac_aas = set("ACDEFGHIKLMNPQRSTVWYX")
        iupac_conform_seqs = self.df["Sequence"].apply(
            lambda seq: set(seq).issubset(iupac_aas)
        )
        if any(~iupac_conform_seqs):
            self.df = self.df.loc[iupac_conform_seqs, :]
            logger.warning(
                f"Sequences are not IUPAC conform. {(~iupac_conform_seqs).sum()} PSMs were dropped."
            )

    def add_protein_ids(self):
        """Add all Protein IDs that matching the sequence.

        Operations are performed inplace on self.df
        """
        peptide_mapper = UPeptideMapper(self.params["database"])
        mapped_peptides = peptide_mapper.map_peptides(self.df["Sequence"].tolist())

        peptide_mappings = [
            merge_and_join_dicts(mapped_peptides[seq], self.DELIMITER)
            for seq in self.df["Sequence"]
        ]

        columns_translations = {
            "start": "Sequence Start",
            "end": "Sequence Stop",
            "post": "Sequence Post AA",
            "id": "Protein ID",
            "pre": "Sequence Pre AA",
        }
        new_columns = pd.DataFrame(peptide_mappings).fillna("")
        new_columns.rename(columns=columns_translations, inplace=True)

        self.df.loc[:, new_columns.columns] = new_columns.values

    def check_enzyme_specificity(self):
        """Check consistency of N/C-terminal cleavage sites.

        Calculates number of missed cleavage sites.
        Operations are performed inplace.
        """
        if self.translated_params["enzyme"]["original_value"] == "nonspecific":
            self.df.loc[:, ["enzN", "enzC"]] = True
            self.df.loc[:, "Missed Cleavages"] = 0
            return None

        enzyme_pattern = self.translated_params["enzyme"]["translated_value"]
        integrity_strictness = self.translated_params[
            "terminal_cleavage_site_integrity"
        ]["translated_value"]

        pren_seq = (
            pd.concat(
                [
                    self.df["Sequence Pre AA"].str.split(rf"{self.DELIMITER}"),
                    self.df["Sequence"].str[:1],
                ],
                axis=1,
            )
            .explode("Sequence Pre AA")
            .sum(axis=1)
        )
        self.df.loc[:, "enzN"] = (
            pren_seq.str.split(fr"{enzyme_pattern}").str[0].str.len() == 1
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
                    self.df["Sequence"].str[-1:],
                    self.df["Sequence Post AA"].str.split("<\\|>"),
                ],
                axis=1,
            )
            .explode("Sequence Post AA")
            .sum(axis=1)
        )
        self.df.loc[:, "enzC"] = (
            postc_seq.str.split(fr"{enzyme_pattern}").str[0].str.len() == 1
        ).groupby(postc_seq.index).agg(integrity_strictness) | (
            postc_seq.str[-1] == "-"
        ).groupby(
            postc_seq.index
        ).agg(
            integrity_strictness
        )

        internal_cuts = self.df["Sequence"].str.split(fr"{enzyme_pattern}")
        self.df.loc[:, "Missed Cleavages"] = (
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
                zip(self.df["Sequence"].values, self.df["Modifications"].values),
                chunksize=1,
            )
        self.df.loc[:, ["uCalc Mass", "Chemical Composition"]] = cc_masses_and_comp
        self.df.loc[:, "uCalc m/z"] = self._calc_mz(
            mass=self.df["uCalc Mass"], charge=self.df["Charge"]
        )
        self.df.loc[:, "Accuracy (ppm)"] = (
            (self.df["Exp m/z"].astype(float) - self.df["uCalc m/z"])
            / self.df["uCalc m/z"]
            * 1e6
        )

    def get_exp_rt_and_mz(self):
        """Experimental mass-to-charge ratios and retention times are added.

        Operations are performed inplace on self.df
        """
        rt_lookup = self._read_rt_lookup_file()
        spec_ids = self.df["Spectrum ID"].astype(int)
        self.df["Retention Time (s)"] = (
            rt_lookup.loc[spec_ids, ["RT", "Unit"]].product(axis=1).to_list()
        )
        self.df["Exp m/z"] = rt_lookup.loc[spec_ids, "Precursor mz"].to_list()

    def add_ranks(self):
        """Ranks are calculated based on the engine scoring column at Spectrum ID level.

        Operations are performed inplace on self.df
        """
        eng_name = self.df["Search Engine"].unique()[0]
        score_col = self.translated_params["validation_score_field"]["translated_value"][
            eng_name
        ]
        top_is_highest = self.translated_params["bigger_scores_better"][
            "translated_value"
        ][eng_name]
        ranking_needs_to_be_ascending = False if top_is_highest is True else True

        self.df.loc[:, "Rank"] = (
            self.df.groupby("Spectrum ID")[score_col]
            .rank(ascending=ranking_needs_to_be_ascending, method="min")
            .astype(int)
        )

    def add_decoy_identity(self):
        """Add boolean decoy state if designated decoy prefix is in Protein IDs.

        Operations are performed inplace on self.df
        """
        decoy_tag = self.params.get("decoy_tag", "decoy_")
        self.df.loc[:, "Is decoy"] = self.df["Protein ID"].str.contains(decoy_tag)

    def sanitize(self):
        """Perform dataframe sanitation steps.

        - Missing raw data locations are replaced by their respective spectrum title identifiers
        - .mgf file extensions are renamed to point to the .mzML files
        - Columns that were not filled in but should exist in the unified format are added and set to None
        - Modifications are sorted alphabetically
        - Columns in the dataframe which could not be properly mapped are removed (warning is raised)
        Operations are performed inplace on self.df
        """
        missing_data_locs = ~(self.df["Raw data location"].str.len() > 0)
        self.df.loc[missing_data_locs, "Raw data location"] = (
            self.df.loc[missing_data_locs, "Spectrum Title"].str.split(".").str[0]
        )
        self.df["Raw data location"] = self.df["Raw data location"].str.replace(
            ".mgf", ".mzML", regex=False
        )

        # Set missing columns to None and reorder columns in standardized manner
        new_cols = self.col_order[~self.col_order.isin(self.df.columns)].to_list()
        self.df.loc[:, new_cols] = None
        self.df = self.df.loc[
            :,
            self.col_order.tolist()
            + sorted(self.df.columns[~self.df.columns.isin(self.col_order)].tolist()),
        ]
        self.df = self.df.astype(self.dtype_mapping)

        # Ensure same order of modifications
        self.df.loc[:, "Modifications"] = (
            self.df["Modifications"]
            .str.split(";")
            .apply(sorted, key=lambda x: x.split(":")[::-1])
            .str.join(";")
        )

        # Ensure there arent any column that should not be
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

    def process_unify_style(self):
        """Combine all additional operations that are needed to calculate new columns and sanitize the dataframe.

        Operations are performed inplace on self.df
        """
        self.df.drop_duplicates(inplace=True, ignore_index=True)
        self.assert_only_iupac_and_missing_aas()
        self.add_protein_ids()
        self.get_exp_rt_and_mz()
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
        self.scan_rt_lookup = self._read_rt_lookup_file()
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
