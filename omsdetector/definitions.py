# -*- coding: utf-8 -*-

# see Zimmermann, N. E. R.; Jain, A.
# Local Structure Order Parameters and Site Fingerprints
# for Quantification of Coordination Environment and Crystal Structure Similarity.
# RSC Adv. 2020, 10 (10), 6063–6081. https://doi.org/10.1039/C9RA07755C.

OP_DEF = {
    4: {
        'names': ['sq_plan', 'sq', 'see_saw_rect', 'tet', 'tri_pyr'],
        'open': [0, 1, 2, 4]
    },
    5: {
        'names': ['pent_plan', 'sq_pyr', 'tri_bipyr'],
        'open': [0, 1]
    },
    6: {
        'names': ['hex_plan_max', 'pent_pyr', 'oct'],
        'open': [0, 1]
    },
    7: {
        'names': ['hex_pyr', 'pent_bipyr'],
        'open': [0]
    },
    8: {
        'names': ['hex_bipyr'],
        'open': None
    }
}
