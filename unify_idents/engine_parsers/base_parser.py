import multiprocessing as mp

import pandas as pd
import uparma
from chemical_composition import ChemicalComposition
from peptide_mapper.mapper import UPeptideMapper
from unimod_mapper.unimod_mapper import UnimodMapper

cc = ChemicalComposition()


def get_cc(seq, mods):
    cc.use(sequence=seq, modifications=mods)
    return cc.mass()


def merge_and_join_dicts(dictlist, delim):
    return {
        key: delim.join([str(d.get(key)) for d in dictlist])
        for key in set().union(*dictlist)
    }


class BaseParser:
    def __init__(self, input_file, params):
        self.input_file = input_file
        if params is None:
            params = {}
        self.params = params
        self.param_mapper = uparma.UParma()

    @classmethod
    def check_parser_compatibility(cls, file):
        return False

    def _read_rt_lookup_file(self):
        rt_lookup = pd.read_csv(self.params["rt_pickle_name"], compression="bz2")
        rt_lookup.set_index("Spectrum ID", inplace=True)
        rt_lookup["Unit"] = rt_lookup["Unit"].replace({"second": 1, "minute": 60})
        return rt_lookup


class __IdentBaseParser(BaseParser):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.DELIMITER = self.params.get("delimiter", "<|>")
        self.PROTON = 1.00727646677
        self.df = None
        self.mod_mapper = UnimodMapper()
        self.params["mapped_mods"] = self.mod_mapper.map_mods(
            mod_list=self.params["modifications"]
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
            "Protein ID": "str",
            "Retention Time (s)": "float32",
            "Exp m/z": "float32",
            "Calc m/z": "float32",
            "uCalc m/z": "float32",
            "uCalc Mass": "float32",
            "Accuracy (ppm)": "float32",
            "Mass Difference": "float32",
            "Sequence Start": "str",
            "Sequence Stop": "str",
            "Sequence Pre AA": "str",
            "Sequence Post AA": "str",
            "Enzyme Specificity": "str",
            "Complies search criteria": "str",
            "Conflicting uparma": "str",
            "Search Engine": "str",
        }
        self.col_order = pd.Series(self.dtype_mapping.keys())

    def _calc_mz(self, mass, charge):
        return (
            mass.astype(float) + (charge.astype(int) * self.PROTON)
        ) / charge.astype(int)

    def _create_mod_dicts(self):
        """Create dict containing meta information about static and variable mods."""
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

    def add_protein_ids(self):
        self.df["Sequence"] = self.df["Sequence"].str.upper()
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
        new_columns = pd.DataFrame(peptide_mappings)
        new_columns.rename(columns=columns_translations, inplace=True)

        self.df.loc[:, new_columns.columns] = new_columns.values

    def calc_masses_and_offsets(self):
        with mp.Pool(self.params.get("cpus", mp.cpu_count() - 1)) as pool:
            cc_masses = pool.starmap(
                get_cc,
                zip(self.df["Sequence"].values, self.df["Modifications"].values),
                chunksize=1,
            )
        self.df.loc[:, "uCalc Mass"] = cc_masses
        self.df.loc[:, "uCalc m/z"] = self._calc_mz(
            mass=self.df["uCalc Mass"], charge=self.df["Charge"]
        )
        self.df.loc[:, "Accuracy (ppm)"] = (
            (self.df["Exp m/z"].astype(float) - self.df["uCalc m/z"])
            / self.df["uCalc m/z"]
            * 1e6
        )

    def get_exp_rt_and_mz(self):
        rt_lookup = self._read_rt_lookup_file()
        spec_ids = self.df["Spectrum ID"].astype(int)
        self.df["Retention Time (s)"] = (
            rt_lookup.loc[spec_ids, ["RT", "Unit"]].product(axis=1).to_list()
        )
        self.df["Exp m/z"] = rt_lookup.loc[spec_ids, "Precursor mz"].to_list()

    def sanitize(self):
        missing_data_locs = ~(self.df["Raw data location"].str.len() > 0)
        self.df.loc[missing_data_locs, "Raw data location"] = (
            self.df.loc[missing_data_locs, "Spectrum Title"].str.split(".").str[0]
        )
        self.df["Raw data location"] = self.df["Raw data location"].str.replace(
            ".mgf", ".mzML", regex=False
        )

        # Set missing columns to None
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
            self.df["Modifications"].str.split(";").apply(sorted).str.join(";")
        )

    def process_unify_style(self):
        self.df.drop_duplicates(inplace=True, ignore_index=True)
        self.add_protein_ids()
        self.calc_masses_and_offsets()
        self.get_exp_rt_and_mz()
        self.sanitize()


class __QuantBaseParser(BaseParser):
    def __init__(self, *args, **kwargs):
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
