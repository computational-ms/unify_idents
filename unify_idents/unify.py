from importlib import import_module
from pathlib import Path
import multiprocessing as mp

import pandas as pd
from chemical_composition import ChemicalComposition
from peptide_mapper.mapper import UPeptideMapper

from unify_idents.engine_parsers.base_parser import BaseParser

cc = ChemicalComposition()


def get_cc(seq, mods):
    cc.use(sequence=seq, modifications=mods)
    return cc.mass()


def merge_and_join_dicts(dictlist, delim):
    return {
        key: delim.join([str(d.get(key)) for d in dictlist])
        for key in set().union(*dictlist)
    }


class Unify:
    """Interface to unify ident outputs from the following engines.

    Args:
        input_file (str): path to file to unify
        params (dict, optional): Description

    Attributes:
        parser (`unify_idents.engine_parsers.base_parser.__BaseParser`): Parser fitting the specified input_file

    """

    def __init__(self, input_file, params={}):
        if not isinstance(input_file, Path):
            self.input_file = Path(input_file)
        else:
            self.input_file = input_file
        self.params = params

        self._parser_classes = []
        self.parser = self._get_parser()
        self.df = None
        self.DELIMITER = self.params.get("delimiter", "<|>")
        self.PROTON = 1.00727646677
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

    def _get_parser(self):
        # Get all files except __init__.pys
        parser_files = (Path(__file__).parent / "engine_parsers").rglob("[!_]*.py")
        # Make paths relative to package
        parser_files = [
            p.relative_to(Path(__file__).parent.parent).as_posix() for p in parser_files
        ]
        # Format and filter for relevant subdirectories
        parser_files = [
            p.replace(".py", "").replace("/", ".")
            for p in parser_files
            if p.count("/") == 3
        ]
        for parser in parser_files:
            import_module(parser)

        for cat in BaseParser.__subclasses__():
            self._parser_classes.extend(cat.__subclasses__())

        for parser in self._parser_classes:
            if parser.check_parser_compatibility(self.input_file) is True:
                return parser(input_file=self.input_file, params=self.params)

        raise IOError(f"No suitable parser found for {self.input_file}.")

    def add_protein_ids(self):
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

        self.df = pd.concat([self.df, new_columns], axis=1)

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
        rt_lookup = pd.read_csv(self.params["rt_pickle_name"], compression="bz2")
        rt_lookup.set_index("Spectrum ID", inplace=True)
        rt_lookup["Unit"] = rt_lookup["Unit"].replace({"second": 1, "minute": 60})
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
        self.df["Sequence"] = self.df["Sequence"].str.upper()

        # Set missing columns to None
        new_cols = self.col_order[~self.col_order.isin(self.df.columns)].to_list()
        self.df.loc[:, new_cols] = None
        self.df = self.df.loc[
            :,
            self.col_order.tolist()
            + sorted(self.df.columns[~self.df.columns.isin(self.col_order)].tolist()),
        ]
        self.df = self.df.astype(self.dtype_mapping)

    def get_dataframe(self):
        self.df = self.parser.unify()
        self.df.drop_duplicates(inplace=True, ignore_index=True)
        self.add_protein_ids()
        self.calc_masses_and_offsets()
        self.get_exp_rt_and_mz()
        self.sanitize()

        return self.df


if __name__ == "__main__":
    mp.freeze_support()
    mp.set_start_method("fork", force=True)
    rt_lookup_path = "/Users/tr341516/PycharmProjects/ursgal2/data/04854_F1_R8_P0109699E13_TMT10_ursgal_lookup.csv.bz2"
    # input_file = "/Users/tr341516/PycharmProjects/ursgal2/data/xtandem_alanine/04854_F1_R8_P0109699E13_TMT10_xtandem_alanine.xml"
    # input_file = "/Users/tr341516/PycharmProjects/ursgal2/data/omssa_2_1_9/04854_F1_R8_P0109699E13_TMT10_omssa_2_1_9.csv"
    input_file = "/Users/tr341516/PycharmProjects/ursgal2/data/msgfplus_2021_03_22/04854_F1_R8_P0109699E13_TMT10_msgfplus_2021_03_22.mzid"
    db_path = "/Users/tr341516/PycharmProjects/ursgal2/data/uniprot_human-ecoli_20180814_IL.fasta"
    obj = Unify(
        input_file,
        params={
            "omssa_mod_dir": Path(__file__).parent.parent / "tests" / "data",
            "rt_pickle_name": rt_lookup_path,
            # "cpus": 1,
            "database": db_path,
            "modifications": [
                {
                    "aa": "M",
                    "type": "opt",
                    "position": "any",
                    "name": "Oxidation",
                },
                {
                    "aa": "C",
                    "type": "fix",
                    "position": "any",
                    "name": "Carbamidomethyl",
                },
                {
                    "aa": "*",
                    "type": "opt",
                    "position": "Prot-N-term",
                    "name": "Acetyl",
                },
            ],
        },
    )
    df = obj.get_dataframe()
    print("done")
