import bz2
import csv
from pathlib import Path

import uparma
from peptide_mapper.mapper import UPeptideMapper
from unimod_mapper.unimod_mapper import UnimodMapper

from unify_idents import UnifiedRow

# get from params
DELIMITER = "<|>"


class __BaseParser:
    def __init__(self, input_file, params=None):
        if params is None:
            self.params = {}
        else:
            self.params = params

        self.param_mapper = uparma.UParma()
        self.peptide_mapper = UPeptideMapper(params["database"])
        self.mod_mapper = UnimodMapper()

        self.map_mods(self.params["Modifications"])

        self.scan_rt_path = self.params.get("scan_rt_lookup_file", None)
        self.scan_rt_lookup = self.read_rt_lookup_file(self.scan_rt_path)

    def file_matches_parser(self):
        # needs to return False to dont be selected as engine parser during `get_parsers`
        return False

    def general_fixes(self, row):
        row["Raw file location"] = row["Spectrum Title"].split(".")[0]
        row["Retention Time (s)"] = float(
            self.scan_rt_lookup[row["Raw file location"]]["scan2rt"][row["Spectrum ID"]]
        )

        # mapping peptides
        starts = []
        ids = []
        stops = []
        pre = []
        post = []
        # we need to convert sequences to uppercase, e.g. omssa reports modified AAs in lowercase
        mapped = self.peptide_mapper.map_peptides([row["Sequence"].upper()])
        for seq, data_list in mapped.items():
            for data in data_list:
                ids.append(data["id"])
                starts.append(str(data["start"]))
                stops.append(str(data["end"]))
                pre.append(str(data["pre"]))
                post.append(str(data["post"]))
        row["Protein ID"] = DELIMITER.join(ids)
        row["Sequence Pre AA"] = DELIMITER.join(pre)
        row["Sequence Post AA"] = DELIMITER.join(post)
        row["Sequence Start"] = DELIMITER.join(starts)
        row["Sequence Stop"] = DELIMITER.join(stops)
        return row

    def read_rt_lookup_file(self, scan_rt_lookup_path):
        with bz2.open(scan_rt_lookup_path, "rt") as fin:
            lookup = {}
            reader = csv.DictReader(fin)
            for line in reader:
                file = Path(line["File"]).stem
                lookup.setdefault(file, {"scan2rt": {}, "rt2scan": {}, "scan2mz": {}})
                scan, rt, mz = (
                    line["Spectrum ID"],
                    line["RT"],
                    line["Precursor mz"],
                )
                lookup[file]["scan2rt"][int(scan)] = float(rt)
                lookup[file]["rt2scan"][float(rt)] = int(scan)
                lookup[file]["scan2mz"][int(scan)] = float(mz)
        return lookup

    def map_mods(self, mods):
        """Maps modifications defined in params["modification"] using unimod.

        Examples:

             >>> [
             ...    "M,opt,any,Oxidation",        # Met oxidation
             ...    "C,fix,any,Carbamidomethyl",  # Carbamidomethylation
             ...    "*,opt,Prot-N-term,Acetyl"    # N-Acteylation
             ... ]

        """
        # TODO remove print functions and replace by logger
        self.params["mods"] = {"fix": [], "opt": []}
        for ursgal_index, mod in enumerate(sorted(mods)):
            mod_params = mod.split(",")
            if len(mod_params) >= 6 or len(mod_params) <= 3:
                print(
                    """
    [ WARNING ] For modifications, please use the ursgal_style:
    [ WARNING ] 'amino_acid,opt/fix,position,Unimod PSI-MS Name'
    [ WARNING ] or
    [ WARNING ] 'amino_acid,opt/fix,position,name,chemical_composition'
    [ WARNING ] Continue without modification {0} """.format(
                        mod
                    )
                )
                print(mod_params)
                continue
            aa = mod_params[0].strip()
            mod_option = mod_params[1].strip()
            pos = mod_params[2].strip()
            unimod = False
            unimod_id = None

            if len(mod_params) == 4:
                try:
                    unimod_id = int(mod_params[3].strip())
                    unimod_name = self.mod_mapper.id2name(unimod_id)
                    mass = self.mod_mapper.id2mass(unimod_id)
                    composition = self.mod_mapper.id2composition(unimod_id)
                    if unimod_name is None:
                        print(
                            """
    [ WARNING ] '{1}' is not a Unimod modification
    [ WARNING ] please change it to a valid Unimod Accession # or PSI-MS Unimod Name
    [ WARNING ] or add the chemical composition hill notation (including 1)
    [ WARNING ] e.g.: H-1N1O2
    [ WARNING ] ursgal_style: 'amino_acid,opt/fix,position,name,chemical_composition'
    [ WARNING ] Continue without modification {0} """.format(
                                mod, unimod_id
                            )
                        )
                        continue
                    unimod = True
                    name = unimod_name
                except:
                    unimod_name = mod_params[3].strip()
                    unimod_id = self.mod_mapper.name2id(unimod_name)
                    mass = self.mod_mapper.name2mass(unimod_name)
                    composition = self.mod_mapper.name2composition(unimod_name)
                    if unimod_id is None:
                        print(
                            """
    [ WARNING ] '{1}' is not a Unimod modification
    [ WARNING ] please change it to a valid PSI-MS Unimod Name or Unimod Accession #
    [ WARNING ] or add the chemical composition hill notation (including 1)
    [ WARNING ] e.g.: H-1N1O2
    [ WARNING ] ursgal_style: 'amino_acid,opt/fix,position,name,chemical_composition'
    [ WARNING ] Continue without modification {0} """.format(
                                mod, unimod_name
                            )
                        )
                        continue
                    unimod = True
                    name = unimod_name

            elif len(mod_params) == 5:
                name = mod_params[3].strip()
                chemical_formula = mod_params[4].strip()
                chemical_composition = ursgal.ChemicalComposition()
                chemical_composition.add_chemical_formula(chemical_formula)
                composition = chemical_composition
                composition_unimod_style = chemical_composition.hill_notation_unimod()
                unimod_name_list = self.mod_mapper.composition2name_list(
                    composition_unimod_style
                )
                unimod_id_list = self.mod_mapper.composition2id_list(
                    composition_unimod_style
                )
                mass = self.mod_mapper.composition2mass(composition_unimod_style)
                for i, unimod_name in enumerate(unimod_name_list):
                    if unimod_name == name:
                        unimod_id = unimod_id_list[i]
                        unimod = True
                        break
                if unimod == False and unimod_name_list != []:
                    print(
                        """
    [ WARNING ] '{0}' is not a Unimod modification
    [ WARNING ] but the chemical composition you specified is included in Unimod.
    [ WARNING ] Please use one of the Unimod names:
    [ WARNING ] {1}
    [ WARNING ] Continue without modification {2} """.format(
                            name, unimod_name_list, mod
                        )
                    )
                    continue
                if unimod == False and unimod_name_list == []:
                    print(
                        """
    [ WARNING ] '{0}' is not a Unimod modification
    [ WARNING ] trying to continue with the chemical composition you specified
    [ WARNING ] This is not working with OMSSA so far""".format(
                            mod,
                        )
                    )
                    mass = chemical_composition._mass()
                    # write new userdefined modifications Xml in unimod style

            mod_dict = {
                "_id": ursgal_index,
                "aa": aa,
                "mass": mass,
                "pos": pos,
                "name": name,
                "composition": composition,
                "org": mod,
                "id": unimod_id,
                "unimod": unimod,
            }
            if mod_dict["unimod"] == False:
                self.mod_mapper.writeXML(mod_dict)

            self.params["mods"][mod_option].append(mod_dict)
        return self.params["mods"]
