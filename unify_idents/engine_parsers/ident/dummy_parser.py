from unify_idents.engine_parsers.base_parser import __IdentBaseParser


class Dummy(__IdentBaseParser):
    def __init__(self, input_file, params=None):
        print("Miss Hoover, I glued my head to my shoulders.")

    def file_matches_parser(self):
        return False
