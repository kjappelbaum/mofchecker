# -*- coding: utf-8 -*-
"""Constants for the MOFChecker"""
# pylint: disable=line-too-long

# see Zimmermann, N. E. R.; Jain, A.
# Local Structure Order Parameters and Site Fingerprints
# for Quantification of Coordination Environment and Crystal Structure Similarity.
# RSC Adv. 2020, 10 (10), 6063–6081. https://doi.org/10.1039/C9RA07755C.

OP_DEF = {
    4: {
        "names": ["sq_plan", "sq", "see_saw_rect", "tet", "tri_pyr"],
        "weights": [0.15, 0.15, 0.1, 0.5, 0.1],
        "open": [0, 1, 2, 4],
    },
    5: {
        "names": ["pent_plan", "sq_pyr", "tri_bipyr"],
        "weights": [1, 1, 0.2],
        "open": [0, 1],
    },
    6: {"names": ["pent_pyr", "oct"], "weights": [0.3, 0.7], "open": [0]},
    7: {
        "names": ["hex_pyr", "pent_bipyr"],
        "weights": [0.7, 0.3],
        "open": [0],
    },
    8: {"names": ["hex_bipyr"], "weights": [1], "open": None},
}

COVALENT_RADII = {
    # Covalent radii revisited -- DOI:10.1039/B801115J
    "H": 0.31,
    "He": 0.28,
    "Li": 1.28,
    "Be": 0.96,
    "B": 0.84,
    "C": 0.76,  # for sp3; sp2 = 0.73; sp = 0.69
    "C_1": 0.69,
    "C_2": 0.73,
    "C_R": 0.73,
    "C_3": 0.76,
    "N": 0.71,
    "O": 0.66,
    "F": 0.57,
    "Ne": 0.58,
    "Na": 1.66,
    "Mg": 1.41,
    "Al": 1.21,
    "Si": 1.11,
    "P": 1.07,
    "S": 1.05,
    "Cl": 1.02,
    "Ar": 1.06,
    "K": 2.03,
    "Ca": 1.76,
    "Sc": 1.7,
    "Ti": 1.6,
    "V": 1.53,
    "Cr": 1.39,
    "Mn": 1.61,  # low spin = 1.39
    "Fe": 1.52,  # low spin = 1.32
    "Co": 1.5,  # low spin = 1.26
    "Ni": 1.24,
    "Cu": 1.32,
    "Zn": 1.22,
    "Ga": 1.22,
    "Ge": 1.2,
    "As": 1.19,
    "Se": 1.2,
    "Br": 1.2,
    "Kr": 1.16,
    "Rb": 2.2,
    "Sr": 1.95,
    "Y": 1.9,
    "Zr": 1.5,  # 1.75  !! temporarily added to correct the UIO cluster bonding problem
    "Nb": 1.64,
    "Mo": 1.54,
    "Tc": 1.47,
    "Ru": 1.46,
    "Rh": 1.42,
    "Pd": 1.39,
    "Ag": 1.45,
    "Cd": 1.44,
    "In": 1.42,
    "Sn": 1.39,
    "Sb": 1.39,
    "Te": 1.38,
    "I": 1.39,
    "Xe": 1.4,
    "Cs": 2.44,
    "Ba": 2.15,
    "La": 2.07,
    "Ce": 2.04,
    "Pr": 2.03,
    "Nd": 2.01,
    "Pm": 1.99,
    "Sm": 1.98,
    "Eu": 1.98,
    "Gd": 1.96,
    "Tb": 1.94,
    "Dy": 1.92,
    "Ho": 1.92,
    "Er": 1.89,
    "Tm": 1.9,
    "Yb": 1.87,
    "Lu": 1.87,
    "Hf": 1.75,
    "Ta": 1.7,
    "W": 1.62,
    "Re": 1.51,
    "Os": 1.44,
    "Ir": 1.41,
    "Pt": 1.36,
    "Au": 1.36,
    "Hg": 1.32,
    "Tl": 1.45,
    "Pb": 1.46,
    "Bi": 1.48,
    "Po": 1.4,
    "At": 1.5,
    "Rn": 1.5,
    "Fr": 2.6,
    "Ra": 2.21,
    "Ac": 2.15,
    "Th": 2.06,
    "Pa": 2,
    "U": 1.96,
    "Np": 1.9,
    "Pu": 1.87,
    "Am": 1.8,
    "Cm": 1.69,
}

EXPECTED_CHECK_VALUES = {
    "has_oms": False,
    "has_carbon": True,
    "has_hydrogen": True,
    "has_atomic_overlaps": False,
    "has_overcoordinated_c": False,
    "has_overcoordinated_n": False,
    "has_overcoordinated_h": False,
    "has_undercoordinated_c": False,
    "has_undercoordinated_n": False,
    "has_metal": True,
    "has_lone_atom": False,
    "has_lone_molecule": False,
    "hash_high_charge": False
}

CHECK_DESCRIPTIONS = {
    "has_oms": "Uses heuristics of order parameter to estimate if there is an uncordinated metal site",
    "has_carbon": "Checks if there is any carbon in the structure",
    "has_hydrogen": "Checks of there is any hydrogen in the structure",
    "has_atomic_overlaps": "Checks if there are atomic overlaps in the structure (estimated based on the adjacency matrix)",
    "has_overcoordinated_c": "Checks if there is any carbon number with coordination number > 4",
    "has_overcoordinated_n": "Checks if there is any nitrogen with coordination number > 4",
    "has_overcoordinated_h": "Checks if there is any carbon with coordination number > 1",
    "has_undercoordinated_c": "Checks with there is any carbon non-linear (i.e., sp2, sp3) carbon with less than two neighbors",
    "has_undercoordinated_n": "Checks if there is a nitrogen that likely misses a hydrogen (e.g., coordinated to a sp2, sp3 carbon)",
    "has_metal": "Checks if there is any metal in the structure",
    "has_lone_atom": "Checks if there is floating atom in the structure",
    "has_lone_molecule": "Checks if there is a floating atom or molecule in the structure",
    "hash_high_charge": "Runs charge equilibration and check if any atom has a high charge (default >3 or <-3)"
}
