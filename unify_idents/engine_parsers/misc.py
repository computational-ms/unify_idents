"""Parser handler."""

import IsoSpecPy as iso
import regex as re
from IsoSpecPy.PeriodicTbl import symbol_to_masses
from chemical_composition import ChemicalComposition
from loguru import logger


def init_custom_cc(function, xml_file_list):
    """Initialize function for multiprocessing by providing 'global' attribute.

    Args:
        function (function): function to be appended with cc attribute
        xml_file_list (list): list of xml files to be passed to ChemicalComposition

    """
    function.cc = ChemicalComposition(unimod_file_list=xml_file_list)


def get_composition_and_mz(seq, mods, charge, exp_mz):
    """Compute hill_notation of any single peptidoform and mass.

    Only the mass of the isotopologue closest to the experimental mass is reported.
    Requires the 'cc' attribute of the function to be set externally.
    Returns None if sequence contains unknown amino acids.
    Args:
        seq (str): peptide sequence
        mods (str): modifications of the peptide sequence, given as "UnimodName:Position"
        charge (int): charge
        exp_mz (float): experimental spectrum mz

    Returns:
        tuple: hill_notation_unimod string, mz

    """
    composition = None
    mz = None
    atom_counts = None
    isotope_masses = None
    isotope_probs = None
    replaced_composition = None
    try:
        get_composition_and_mz.cc.use(sequence=seq, modifications=mods)
        composition = get_composition_and_mz.cc.hill_notation_unimod()
        static_isotopes = re.findall(r"(?<=\))(\d+)(\w+)(?:\()(\d+)", composition)
        if len(static_isotopes) != 0:
            atom_counts = []
            isotope_masses = []
            for isotope in static_isotopes:
                mass, element, number = isotope
                replaced_composition = composition.replace(
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
        isotopologue_mzs = list(
            iso.IsoThreshold(
                formula=formula,
                threshold=0.001,
                charge=charge,
                get_confs=True,
                atomCounts=atom_counts,
                isotopeMasses=isotope_masses,
                isotopeProbabilities=isotope_probs,
            ).masses
        )
        # Report only most accurate mass
        mz = min(
            isotopologue_mzs,
            key=lambda x: abs(exp_mz - x),
        )
    except (KeyError, Exception):
        logger.warning(f"Could not calculate mz for {seq}#{mods}")
        pass
    return composition, mz