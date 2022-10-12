"""Parser handler."""

import uparma


class BaseParser:
    """Base class of all parser types."""

    def __init__(self, input_file, params, immutable_peptides=None):
        """Initialize parser.

        Args:
            input_file (str): path to input file
            params (dict): ursgal param dict
            immutable_peptides (list, optional): list of immutable peptides
        """
        self.input_file = input_file
        self.immutable_peptides = immutable_peptides
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
