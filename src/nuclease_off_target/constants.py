# -*- coding: utf-8 -*-
"""Genomic sequences."""
from typing import Dict
from typing import Union


SECONDS_BETWEEN_UCSC_REQUESTS = 10
# for displaying vertical alignments
VERTICAL_ALIGNMENT_MATCH_CHARACTER = " "
VERTICAL_ALIGNMENT_MISMATCH_CHARACTER = "X"
ALIGNMENT_GAP_CHARACTER = "-"
VERTICAL_ALIGNMENT_DNA_BULGE_CHARACTER = "+"
VERTICAL_ALIGNMENT_RNA_BULGE_CHARACTER = "-"
SEPARATION_BETWEEN_GUIDE_AND_PAM = " "

CAS_VARIETIES: Dict[str, Dict[str, Union[int, str, Dict[int, Union[int, float]]]]] = {
    "Sa": {
        "PAM": "NNGRRT",
        "cut_site_relative_to_pam": -3,
        "mismatch-penalties-starting-from-PAM": {
            0: 6,
            1: 5,
            2: 4,
            3: 3,
            21: 0.1,
            20: 0.1,
            19: 0.12,
            18: 0.13,
            17: 0.15,
            16: 0.17,
            15: 0.19,
            14: 0.21,
            13: 0.23,
            12: 0.27,
            11: 0.35,
            10: 0.5,
            9: 0.7,
            8: 0.8,
            7: 1.1,
            6: 1.3,
            5: 1.9,
            4: 2.3,
        },
    }
}
