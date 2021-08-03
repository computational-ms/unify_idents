from unify_idents.engine_parsers.base_parser import __QuantBaseParser


class Dummy(__QuantBaseParser):
    def __init__(self, input_file, params=None):
        print("Test")

    def file_matches_parser(self):
        return False
