"""Parser handler."""

import IsoSpecPy as iso
import numpy as np
import regex as re
from IsoSpecPy.PeriodicTbl import symbol_to_masses
from chemical_composition import ChemicalComposition
from loguru import logger


def trunc(values, decs=0):
    """Truncate  float to `number of decimals.

    Args:
        decs (int): Truncating precision
    """
    return np.trunc(np.around(values, decimals=decs) * 10**decs) / (10**decs)

    # return np.trunc(values * 10**decs) / (10**decs)


def init_custom_cc(function, xml_file_list, proton):
    """Initialize function for multiprocessing by providing 'global' attribute.

    Args:
        function (function): function to be appended with cc attribute
        xml_file_list (list): list of xml files to be passed to ChemicalComposition
        proton (float): proton mass
    """

    def _calc_mz(mass, charge):
        return (mass + (charge * proton)) / charge

    function.cc = ChemicalComposition(unimod_file_list=xml_file_list)
    function.calc_mz = _calc_mz


def get_composition_and_mass_and_accuracy(seq, mods, charge, exp_mz):
    """Compute hill_notation of any single peptidoform, mass, and accuracy.

    Only the accuracy of the isotopologue closest to the experimental mass is reported.
    Requires the 'cc' attribute of the function to be set externally.
    Returns None if sequence contains unknown amino acids.

    Args:
        seq (str): peptide sequence
        mods (str): modifications of the peptide sequence, given as "UnimodName:Position"
        charge (int): charge
        exp_mz (float): experimental spectrum mz
    Returns:
        tuple: hill_notation_unimod string, mass, accuracy
    """
    composition = None
    c12_mass = None
    isotopologue_acc = None
    atom_counts = None
    isotope_masses = None
    isotope_probs = None
    try:
        get_composition_and_mass_and_accuracy.cc.use(sequence=seq, modifications=mods)
        composition = get_composition_and_mass_and_accuracy.cc.hill_notation_unimod()
        replaced_composition = composition
        c12_mass = get_composition_and_mass_and_accuracy.cc.mass()
        static_isotopes = re.findall(r"(?<=\))(\d+)(\w+)(?:\()(\d+)", composition)
        if len(static_isotopes) != 0:
            atom_counts = []
            isotope_masses = []
            for isotope in static_isotopes:
                mass, element, number = isotope
                replaced_composition = replaced_composition.replace(
                    f"{mass}{element}({number})", ""
                )
                isotope_masses.append(
                    [
                        [
                            m
                            for m in symbol_to_masses[element]
                            if str(round(m)).startswith(mass)
                        ][0]
                    ]
                )
                atom_counts.append(int(number))
            isotope_probs = len(atom_counts) * [[1.0]]
        if replaced_composition is not None:
            formula = replaced_composition
        else:
            formula = composition
        formula = formula.replace("(", "").replace(")", "")
        isotopologue_masses = list(
            iso.IsoThreshold(
                formula=formula,
                threshold=0.02,
                charge=1,
                get_confs=True,
                atomCounts=atom_counts,
                isotopeMasses=isotope_masses,
                isotopeProbabilities=isotope_probs,
            ).masses
        )
        isotopologue_mzs = [
            get_composition_and_mass_and_accuracy.calc_mz(mass, charge)
            for mass in isotopologue_masses
        ]
        # Report only most accurate mass
        isotopologue_mz = min(
            isotopologue_mzs,
            key=lambda x: abs(exp_mz - x),
        )
        isotopologue_acc = (exp_mz - isotopologue_mz) / isotopologue_mz * 1e6
    except (KeyError, Exception):
        logger.warning(f"Could not calculate mz for {seq}#{mods}")
        pass
    return composition, c12_mass, isotopologue_acc
