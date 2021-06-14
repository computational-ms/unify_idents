import csv
from unify_idents import (
    UnifiedRow,
)
import uparma


"""All engines
        * Retention Time (s) is correctly set using _ursgal_lookup.pkl
          During mzML conversion to mgf the retention time for every spec
          is stored in a internal lookup and used later for setting the RT.
        * All modifications are checked if they were given in
          params['modifications'], converted to the name that was given
          there and sorted according to their position.
        * Fixed modifications are added in 'Modifications', if not reported
          by the engine.
        * The monoisotopic m/z for for each line is calculated (uCalc m/z),
          since not all engines report the monoisotopic m/z
        * Mass accuracy calculation (in ppm), also taking into account that
          not always the monoisotopic peak is picked
        * Rows describing the same PSM (i.e. when two proteins share the
          same peptide) are merged to one row.
"""


class __BaseParser:
    def __init__(self, input_file, params=None):
        self.param_mapper = uparma.UParma()
        self.peptide_mapper = None
        # should it also have scan_rt_lookup?

    def file_matches_parser(self):
        # needs to return False to dont be selected as engine parser during `get_parsers`
        return False

    def general_fixes(self, row):
        # format mods
        # peptide mapping
        return row
