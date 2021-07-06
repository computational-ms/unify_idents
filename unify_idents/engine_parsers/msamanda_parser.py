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
        self.column_mapping = self.get_column_names()

    
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
    @classmethod
    def file_matches_parser(cls, file):
        # TODO implement file sensing
        # use get column names
        msamanda_version = "#version: 2.0.0.17442"
        with open(file) as fh:
            reader = csv.DictReader(fh)
            if reader.fieldnames[0] == msamanda_version:
                ret_val = True
            else:
                ret_val = False   
        return ret_val

    def __iter__(self):
        return self
    
    def __next__(self):
        n = next(self.reader)
        u = self._unify_row(n)
        return u
    
    # def _unify_row(self, row):
    
    #     new_row = {}
        # for unify_name, omssa_name in self.column_mapping.items():
        #     new_row[unify_name] = row[omssa_name]
        # for col in self.cols_to_remove:
        #     del new_row[col]
        # for col in self.cols_to_add:
        #     new_row[col] = ""
        # new_row["Spectrum ID"] = int(new_row["Spectrum Title"].split(".")[1])
        # new_row["Search Engine"] = "omssa_2_1_9"
    
        # modstring = self.create_mod_string(new_row)
        # new_row["Modifications"] = modstring
        # new_row = self.general_fixes(new_row)
    
        # return UnifiedRow(**new_row)
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
    def get_column_names(self):
        # create own uparma mapper
        headers = self.param_mapper.get_default_params(style=self.style)[
            "header_translations"
        ]["translated_value"]
        return headers
    

# def output_cleanup(udict):
#     cached_msamanada_output = []
#     result_file = open(f"{udict['output_filenames'][0]}", "r")
#     csv_dict_reader_object = csv.DictReader(
#         (row for row in result_file if not row.startswith("#")), delimiter="\t"
#     )
#     headers = csv_dict_reader_object.fieldnames
#     header_translations = udict["translated_cparameters"]["header_translations"][
#         "translated_value"
#     ]
#     header_translations = {value: key for key, value in header_translations.items()}
#     translated_headers = []
#     for header in headers:
#         translated_headers.append(header_translations.get(header, header))
#     translated_headers += ["Raw data location"]
#     logger.info("[ PARSING  ] Loading unformatted MS Amanda results ...")
#     for line_dict in csv_dict_reader_object:
#         cached_msamanada_output.append(line_dict)
#     logger.info("Loading unformatted MS Amanda results done!")
#     result_file.close()

#     # self.params["output_file"] = self.params["output_file"].replace("tsv", "csv")
#     if sys.platform == "win32":
#         lineterminator = "\n"
#     else:
#         lineterminator = "\r\n"
#     with open(f"{udict['output_filenames'][0]}", "w") as result_file:
#         csv_dict_writer_object = csv.DictWriter(
#             result_file,
#             fieldnames=translated_headers,
#             lineterminator=lineterminator,
#         )
#         csv_dict_writer_object.writeheader()
#         logger.info("Writing MS Amanda results, this can take a while...")
#         csv_write_list = []
#         total_docs = len(cached_msamanada_output)
#         for cache_pos, m in enumerate(cached_msamanada_output):
#             tmp = {}
#             for header in headers:
#                 translated_header = header_translations.get(header, header)
#                 tmp[translated_header] = m[header]
#             tmp["Sequence"] = tmp["Sequence"].upper()

#             if cache_pos % 500 == 0:
#                 logger.info(
#                     "[ INFO ] Processing line number:    {0}/{1}".format(
#                         cache_pos, total_docs
#                     ),
#                     end="\r",
#                 )

#             dict_2_write = copy.deepcopy(tmp)
#             translated_mods = []
#             # N-Term(Acetyl|42.010565|fixed);M1(Oxidation|15.994915|fixed);M23(Oxidation|15.994915|fixed)
#             if dict_2_write["Modifications"] != "":
#                 splitted_Modifications = dict_2_write["Modifications"].split(";")
#                 for mod in splitted_Modifications:

#                     (
#                         position_or_aa_and_pos_unimod_name,
#                         mod_mass,
#                         fixed_or_opt,
#                     ) = mod.split("|")
#                     (
#                         position_or_aa_and_pos,
#                         unimod_name,
#                     ) = position_or_aa_and_pos_unimod_name.split("(")
#                     position_or_aa_and_pos = position_or_aa_and_pos.strip()
#                     unimod_name = unimod_name.strip()

#                     if position_or_aa_and_pos.upper() == "N-TERM":
#                         position = 0
#                     else:
#                         position = position_or_aa_and_pos[1:]

#                     translated_mods.append("{0}:{1}".format(unimod_name, position))

#             dict_2_write["Modifications"] = ";".join(translated_mods)
#             dict_2_write["Raw data location"] = udict["ufiles"][0]

#             # protein_id = tmp['proteinacc_start_stop_pre_post_;']

#             csv_write_list.append(dict_2_write)
#         duplicity_buffer = set()
#         for final_dict_2_write in csv_write_list:
#             duplicity_key = (
#                 final_dict_2_write["Sequence"],
#                 final_dict_2_write["Modifications"],
#                 final_dict_2_write["proteinacc_start_stop_pre_post_;"],
#                 final_dict_2_write["Spectrum ID"],
#             )
#             if duplicity_key not in duplicity_buffer:
#                 csv_dict_writer_object.writerow(final_dict_2_write)
#                 duplicity_buffer.add(duplicity_key)
#     logger.info("Writing MS Amanda results done!")
#     return

