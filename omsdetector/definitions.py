# -*- coding: utf-8 -*-

# see Zimmermann, N. E. R.; Jain, A.
# Local Structure Order Parameters and Site Fingerprints
# for Quantification of Coordination Environment and Crystal Structure Similarity.
# RSC Adv. 2020, 10 (10), 6063–6081. https://doi.org/10.1039/C9RA07755C.

OP_DEF = {
    4: {
        'names': ['sq_plan', 'sq', 'see_saw_rect', 'tet', 'tri_pyr'],
        'weights': [0.1, 0.1, 0.1, 0.6, 0.1],
        'open': [0, 1, 2, 4]
    },
    5: {
        'names': ['pent_plan', 'sq_pyr', 'tri_bipyr'],
        'weights': [1, 1, 1],
        'open': [0, 1]
    },
    6: {
        'names': ['hex_plan_max', 'pent_pyr', 'oct'],
        'weights': [0.1, 0.1, 0.8],
        'open': [0, 1]
    },
    7: {
        'names': ['hex_pyr', 'pent_bipyr'],
        'weights': [0.7, 0.3],
        'open': [0],
    },
    8: {
        'names': ['hex_bipyr'],
        'weights': [1],
        'open': None
    }
}
