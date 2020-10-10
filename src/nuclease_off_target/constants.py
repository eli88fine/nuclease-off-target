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

CAS_VARIETIES: Dict[str, Dict[str, Union[int, str]]] = {
    "Sa": {"PAM": "NNGRRT", "cut_site_relative_to_pam": -3}
}
