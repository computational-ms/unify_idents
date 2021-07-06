#!/usr/bin/env python
import csv

import uparma

from unify_idents import UnifiedRow
from unify_idents.engine_parsers.base_parser import __BaseParser

import xml.etree.ElementTree
from pathlib import Path
from loguru import logger


class MSamandaParser(__BaseParser):
    def __init__(self, input_file, params=None):
        super().__init__(input_file, params)
        if params is None:
            params = {}
        self.params = params
        self.input_file = input_file

        try:
            self.reader = csv.DictReader(open(input_file))
        except:
            self.reader = None

        self.style = "msamanda_style_1"
        # self.column_mapping = self.get_column_names()

    #     # self.omssa_mod_corrections = {
    #     #     # this dict holds corrections of wrong OMSSA to unimod assignments
    #     #     "TMT 6-plex on K": {
    #     #         "unimod_id": "737",
    #     #         "ommsa_unimod_id": "738",  # this is TMT duplex in unimod
    #     #         "unimod_name": "TMT6plex",
    #     #     },
    #     #     "TMT 6-plex on n-term peptide": {
    #     #         "unimod_id": "737",
    #     #         "ommsa_unimod_id": "738",  # this is TMT duplex in unimod
    #     #         "unimod_name": "TMT6plex",
    #     #         "aa_targets": ["N-term"],  # override 'X' in OMSSA mods xml
    #     #     },
    #     #     "TMT duplex on K": {
    #     #         "unimod_id": "738",
    #     #         "ommsa_unimod_id": "738",  # this is TMT duplex in unimod
    #     #         "unimod_name": "TMT2plex",
    #     #     },
    #     #     "TMT duplex on n-term peptide": {
    #     #         "unimod_id": "738",
    #     #         "ommsa_unimod_id": "738",  # this is TMT duplex in unimod
    #     #         "unimod_name": "TMT2plex",
    #     #         "aa_targets": ["N-term"],  # override 'X' in OMSSA mods xml
    #     #     },
    #     # }
    #
    #     self.cols_to_remove = [
    #         "proteinacc_start_stop_pre_post_;",
    #         "Start",
    #         "Stop",
    #         "NIST score",
    #         "gi",
    #         "Accession",
    #     ]
    #
    #     self.cols_to_add = [
    #         "uCalc m/z",
    #         "uCalc Mass",
    #         "Retention Time (s)",
    #         "Accuracy (ppm)",
    #         "Mass Difference",
    #         "Protein ID",
    #         "Sequence Start",
    #         "Sequence Stop",
    #         "Sequence Pre AA",
    #         "Sequence Post AA",
    #         "Enzyme Specificity",
    #         "Complies search criteria",
    #         "Conflicting uparam",
    #         "Search Engine",
    #     ]
    #     if self.reader is not None:
    #         self.create_mod_lookup()
    #
    # @classmethod
    # def file_matches_parser(cls, file):
    #     # TO.DO implement file sensing
    #     # use get column names
    #     fn = [
    #         "Spectrum number",
    #         " Filename/id",
    #         " Peptide",
    #         " E-value",
    #         " Mass",
    #         " gi",
    #         " Accession",
    #         " Start",
    #         " Stop",
    #         " Defline",
    #         " Mods",
    #         " Charge",
    #         " Theo Mass",
    #         " P-value",
    #         " NIST score",
    #     ]
    #     field_set = set([a.strip() for a in fn])
    #     with open(file) as fh:
    #         reader = csv.DictReader(fh)
    #         if set([f.strip() for f in reader.fieldnames]) == field_set:
    #             ret_val = True
    #         else:
    #             ret_val = False
    #     return ret_val
    #
    # def __iter__(self):
    #     return self
    #
    # def __next__(self):
    #     n = next(self.reader)
    #     u = self._unify_row(n)
    #     return u
    #
    # def _unify_row(self, row):
    #
    #     new_row = {}
    #     for unify_name, omssa_name in self.column_mapping.items():
    #         new_row[unify_name] = row[omssa_name]
    #     for col in self.cols_to_remove:
    #         del new_row[col]
    #     for col in self.cols_to_add:
    #         new_row[col] = ""
    #     new_row["Spectrum ID"] = int(new_row["Spectrum Title"].split(".")[1])
    #     new_row["Search Engine"] = "omssa_2_1_9"
    #
    #     modstring = self.create_mod_string(new_row)
    #     new_row["Modifications"] = modstring
    #     new_row = self.general_fixes(new_row)
    #
    #     return UnifiedRow(**new_row)
    #
    # def create_mod_string(self, new_row):
    #     fixed_mods = []
    #     for pos, aa in enumerate(new_row["Sequence"]):
    #         pos += 1  # start counting at 1, 0 is N-term
    #         for mod in self.params["mods"]["fix"]:
    #             if aa == mod["aa"]:
    #                 fixed_mods.append(f"{mod['name']}:{pos}")
    #     f_mod_string = ";".join(fixed_mods)
    #
    #     translated_mods = []
    #     if new_row["Modifications"] != "":
    #         splitted_Modifications = new_row["Modifications"].split(",")
    #         for mod in splitted_Modifications:
    #             omssa_name, position = mod.split(":")
    #             omssa_name = omssa_name.strip()
    #             position = position.strip()
    #             unimod_name = self.lookups[omssa_name]["name"]
    #             if position.strip() == "1":
    #                 # print( self.lookups[ omssa_name ] )
    #                 for target in self.lookups[omssa_name]["aa_targets"]:
    #                     if "N-TERM" in target.upper():
    #                         position = "0"
    #             translated_mods.append("{0}:{1}".format(unimod_name, position))
    #
    #     # join fixed and variable mods
    #
    #     fix_mods = [mod.split(":") for mod in f_mod_string.split(";") if mod != ""]
    #     opt_mods = [mod.split(":") for mod in translated_mods if mod != ""]
    #     all_mods = sorted(fix_mods + opt_mods, key=lambda x: float(x[1]))
    #     all_mods = ";".join([":".join(m) for m in all_mods])
    #     return all_mods
    #
    # def get_column_names(self):
    #     # create own uparma mapper
    #     headers = self.param_mapper.get_default_params(style=self.style)[
    #         "header_translations"
    #     ]["translated_value"]
    #     return headers
    #
    # def _load_omssa_xml(self):
    #     """Parsing through omssa mods to map omssa mods on unimods"""
    #     self.omssa_mod_mapper = {}
    #
    #     def _create_empty_tmp():
    #         tmp = {
    #             "aa_targets": [],
    #         }
    #         return tmp
    #
    #     tmp = _create_empty_tmp()
    #     xml_path = Path(self.params["omssa_mod_dir"])
    #
    #     xml_names = ["mods.xml", "usermods.xml"]
    #     for xml_name in xml_names:
    #         omssa_xml = Path(xml_path) / xml_name
    #         # self.print_info(
    #         #     "Parsing omssa xml ({0})".format(omssa_xml), caller="__ini__"
    #         # )
    #         for event, element in xml.etree.ElementTree.iterparse(omssa_xml):
    #             if element.tag.endswith("MSModSpec_residues_E"):
    #                 tmp["aa_targets"].append(element.text)
    #
    #             elif element.tag.endswith("MSMod"):
    #                 tmp["omssa_id"] = element.text
    #                 # tmp['MSMod'] = element.text # OMSSA ID!
    #             elif element.tag.endswith("MSModSpec_psi-ms"):
    #                 tmp["unimod_name"] = element.text
    #             elif element.tag.endswith("MSModSpec_unimod"):
    #                 tmp["unimod_id"] = element.text
    #                 # tmp['MSModSpec_psi-ms'] = element.text # UNIMOD Name
    #             elif element.tag.endswith("MSModSpec_name"):
    #                 tmp["omssa_name"] = element.text
    #                 additional = []
    #                 if "protein" in tmp["omssa_name"]:
    #                     additional.append("Prot")
    #                 if "n-term" in tmp["omssa_name"]:
    #                     additional.append("N-term")
    #                 elif "c-term" in tmp["omssa_name"]:
    #                     additional.append("C-term")
    #                 if len(additional) > 0:
    #                     tmp["aa_targets"].append("-".join(additional))
    #
    #             elif element.tag.endswith("MSModSpec"):
    #                 lookup_field = "unimod_id"
    #                 try:
    #                     l_value = tmp[lookup_field]
    #                 except:
    #                     tmp["aa_targets"] = []
    #                     continue
    #                 if tmp["omssa_name"] in self.omssa_mod_corrections.keys():
    #                     l_value = self.omssa_mod_corrections[tmp["omssa_name"]][
    #                         "unimod_id"
    #                     ]
    #                     # for TMT mods OMSSA writes an 'X' as amino acid target, this breaks later code...
    #                     if (
    #                         "aa_targets"
    #                         in self.omssa_mod_corrections[tmp["omssa_name"]].keys()
    #                     ):
    #                         tmp["aa_targets"] = self.omssa_mod_corrections[
    #                             tmp["omssa_name"]
    #                         ]["aa_targets"]
    #                 if l_value not in self.omssa_mod_mapper.keys():
    #                     self.omssa_mod_mapper[l_value] = {}
    #                 self.omssa_mod_mapper[l_value][tmp["omssa_id"]] = {
    #                     "aa_targets": tmp["aa_targets"],
    #                     "omssa_name": tmp["omssa_name"],
    #                 }
    #
    #                 tmp = _create_empty_tmp()
    #     return
    #
    # def create_mod_lookup(self):
    #     self._load_omssa_xml()
    #     self.lookups = {}
    #     for param_key in ["_fixed_mods", "_opt_mods"]:
    #         mod_type = param_key[1:4]
    #         modifications = ""
    #         self.params[param_key] = ""
    #         for mod in self.params["mods"][mod_type]:
    #             unimod_id_does_not_exist = False
    #             aa_can_not_be_mapped = True
    #             if mod["id"] not in self.omssa_mod_mapper.keys():
    #                 unimod_id_does_not_exist = True
    #             else:
    #                 if mod["aa"] == "*":
    #                     search_target = [
    #                         mod["pos"],
    #                     ]
    #                 else:
    #                     search_target = [
    #                         mod["aa"],
    #                     ]
    #                 for omssa_id in self.omssa_mod_mapper[mod["id"]].keys():
    #                     if (
    #                         search_target
    #                         == self.omssa_mod_mapper[mod["id"]][omssa_id]["aa_targets"]
    #                     ):
    #                         modifications += "{0},".format(omssa_id)
    #                         aa_can_not_be_mapped = False
    #                         omssa_name = self.omssa_mod_mapper[mod["id"]][omssa_id][
    #                             "omssa_name"
    #                         ]
    #                         self.lookups[omssa_name] = {
    #                             "name": mod["name"],
    #                             "aa_targets": self.omssa_mod_mapper[mod["id"]][
    #                                 omssa_id
    #                             ]["aa_targets"],
    #                             "omssa_id": omssa_id,
    #                             "id": mod["id"],
    #                         }
    #             if unimod_id_does_not_exist or aa_can_not_be_mapped:
    #                 logger.warning(
    #                     """
    # The combination of modification name and aminoacid is not supported by
    # OMSSA. Continuing without modification: {0}
    #                 """.format(
    #                         mod
    #                     ),
    #                 )
    #                 continue
    #
    #         self.params[param_key] = modifications.strip(",")
