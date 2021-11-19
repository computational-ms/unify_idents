from importlib import import_module
from pathlib import Path

from unify_idents.engine_parsers.base_parser import BaseParser


class Unify:
    """Interface to unify ident outputs from different engines.

    Args:
        input_file (str): path to file to unify
        params (dict, optional): Description

    Attributes:
        parser (`unify_idents.engine_parsers.base_parser.__BaseParser`): Parser fitting the specified input_file

    """

    def __init__(self, input_file, params={}):
        """
        Args:
            input_file (str): path to input file
            params (dict): ursgal param dict
        """
        if not isinstance(input_file, Path):
            self.input_file = Path(input_file)
        else:
            self.input_file = input_file
        self.params = params

        self._parser_classes = []
        self.parser = self._get_parser()
        self.df = None
        self.PROTON = 1.00727646677

    def _get_parser(self):
        """
        Checks input file / parser compatibility and init matching parser in self.parser.
        Raises error if no matching parser can be found.
        """
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

    def get_dataframe(self):
        """
        Computes and returns a unified dataframe using the input file with specified parameters.
        Returns:
            self.df (pd.DataFrame): unified dataframe

        """
        self.df = self.parser.unify()

        return self.df


if __name__ == "__main__":
    rt_lookup_path = "/Users/tr341516/PycharmProjects/ursgal2/data/04854_F1_R8_P0109699E13_TMT10_ursgal_lookup.csv.bz2"
    input_file = "/Users/tr341516/PycharmProjects/ursgal2/data/xtandem_vengeance_8bc71feada0641dac433424377a4b36b/04854_F1_R8_P0109699E13_TMT10_xtandem_vengeance.xml"
    # input_file = "/Users/tr341516/PycharmProjects/ursgal2/data/msgfplus_2021_03_22_53574949adc37d787f5ae9934c14cf36/04854_F1_R8_P0109699E13_TMT10_msgfplus_2021_03_22.mzid"
    # input_file = "/Users/tr341516/PycharmProjects/ursgal2/data/msgfplus_2021_03_22/04854_F1_R8_P0109699E13_TMT10_msgfplus_2021_03_22.mzid"
    # input_file = "/Users/tr341516/PycharmProjects/ursgal2/data/msfragger_3/04854_F1_R8_P0109699E13_TMT10_msfragger_3.tsv"
    # input_file = "/Users/tr341516/PycharmProjects/unify_idents/tests/data/test_Creinhardtii_QE_pH11_msamanda_2_0_0_17442.csv"
    # input_file = "/Users/tr341516/PycharmProjects/ursgal2/data/comet_2020_01_4/04854_F1_R8_P0109699E13_TMT10_comet_2020_01_4.mzid"
    # input_file = "/Users/tr341516/Downloads/07634_F1_R3_P0215547J46_TMT11_thermo_raw_file_parser_1_1_11_comet_2020_01_4.mzid"
    # input_file = "/Users/tr341516/PycharmProjects/unify_idents/tests/data/flash_lfq_1_2_0_quantified_peaks.tsv"
    db_path = "/Users/tr341516/PycharmProjects/ursgal2/data/uniprot_human-ecoli_20180814_IL.fasta"
    obj = Unify(
        input_file,
        params={
            "omssa_mod_dir": Path(__file__).parent.parent / "tests" / "data",
            "rt_pickle_name": rt_lookup_path,
            "Raw data location": "/Users/tr341516/PycharmProjects/ursgal2/data/04854_F1_R8_P0109699E13_TMT10.mzML",
            "database": db_path,
            "modifications": [
                {
                    "aa": "M",
                    "type": "opt",
                    "position": "any",
                    "name": "Oxidation",
                },
                {
                    "aa": "*",
                    "type": "opt",
                    "position": "Prot-N-term",
                    "name": "Acetyl",
                },
                {
                    "aa": "*",
                    "type": "opt",
                    "position": "N-term",
                    "name": "TMT6plex",
                },
                {
                    "aa": "C",
                    "type": "fix",
                    "position": "any",
                    "name": "Carbamidomethyl",
                },
                {
                    "aa": "K",
                    "type": "fix",
                    "position": "any",
                    "name": "TMT6plex",
                },
            ],
        },
    )
    df = obj.get_dataframe()
    print("done")
