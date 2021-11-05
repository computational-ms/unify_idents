import pandas as pd
import uparma
from unimod_mapper.unimod_mapper import UnimodMapper


class BaseParser:
    def __init__(self, input_file, params):
        self.PROTON = 1.00727646677
        self.input_file = input_file
        if params is None:
            params = {}
        self.params = params
        self.param_mapper = uparma.UParma()
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

    @classmethod
    def check_parser_compatibility(cls, file):
        return False

    def _calc_mz(self, mass, charge):
        return (
            mass.astype(float) + (charge.astype(int) * self.PROTON)
        ) / charge.astype(int)


class __IdentBaseParser(BaseParser):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.mod_mapper = UnimodMapper()
        self.params["mapped_mods"] = self.mod_mapper.map_mods(
            mod_list=self.params["modifications"]
        )
        self.mod_dict = self._create_mod_dicts()

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
