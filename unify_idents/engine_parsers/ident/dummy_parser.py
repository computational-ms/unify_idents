"""Dummy parser."""
from unify_idents.engine_parsers.base_parser import IdentBaseParser


class Dummy(IdentBaseParser):
    """File parser for Dummy."""

    def __init__(self, *args, **kwargs):
        """Initialize parser.

        Reads in data file and provides mappings.
        """
        super().__init__(*args, **kwargs)
        print("Miss Hoover, I glued my head to my shoulders.")

    @classmethod
    def check_parser_compatibility(cls, file):
        """Assert compatibility between file and parser.

        Args:
            file (str): path to input file

        Returns:
            bool: True if parser and file are compatible

        """
        return False

    def unify(self):
        """Primary method to read and unify engine output.

        Returns:
            self.df (pd.DataFrame): unified dataframe
        """
        return None
