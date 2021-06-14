#!/usr/bin/env python3


class UnifiedRow:
    def __init__(self, **kwargs):
        self.columns = kwargs.keys()
        self.data = kwargs
        self.col_order = [
            "Spectrum Title",
            "Spectrum ID",
            "Sequence",
            "Modifications",
            "Charge",
            "Protein ID",
            "Retention Time (s)",
            "Exp m/z",
            "Calc m/z",
            "uCalc_m/z",
            "uCalc Mass",
            "Accuracy (ppm)",
            "Mass Difference",
            "Sequence Start",
            "Sequence Stop",
            "Sequence Pre AA",
            "Sequence Post AA",
            "Enzyme Specificity",
            "Complies search criteria",
            "Conflicting uparma",
            "Search Engine",
        ]
        self.string_repr = None

    def __str__(self):
        if self.string_repr is None:
            self.string_repr = []
            for col in self.col_order:
                self.string_repr.append(self.data.get(col, ""))
            self.string_repr = ", ".join(self.string_repr)
        return self.string_repr

    def __repr__(self):
        return self.__str__()

    def calc_mz(self):
        # use chemical composition
        # using pyqms would require each row to initialize a isotopologue lib object ...
        pass

    def cast_types(self):
        # TODO cast to correct types
        pass
