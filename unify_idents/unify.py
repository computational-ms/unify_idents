"""Unify handler."""
from importlib import import_module
from pathlib import Path

from unify_idents.engine_parsers.base_parser import BaseParser


class Unify:
    """Interface to unify ident outputs from different engines.

    Args:
        input_file (str): path to file to unify
        params (dict): Description

    Attributes:
        parser (`unify_idents.engine_parsers.base_parser.__BaseParser`): Parser fitting the specified input_file

    """

    def __init__(self, input_file, params):
        """Initialize unifier.

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
        """Check input file / parser compatibility and init matching parser in self.parser.

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
        """Compute and returns a unified dataframe using the input file with specified parameters.

        Returns:
            self.df (pd.DataFrame): unified dataframe

        """
        self.df = self.parser.unify()

        return self.df
