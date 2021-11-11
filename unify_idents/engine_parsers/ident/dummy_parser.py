from unify_idents.engine_parsers.base_parser import __IdentBaseParser


class Dummy(__IdentBaseParser):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        print("Miss Hoover, I glued my head to my shoulders.")

    @classmethod
    def check_parser_compatibility(cls, file):
        return False

    def unify(self):
        return None
