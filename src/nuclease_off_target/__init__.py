# -*- coding: utf-8 -*-
"""Docstring."""
from . import genomic_sequence
from .constants import ALIGNMENT_GAP_CHARACTER
from .constants import ALIGNMENT_MATCH_CHARACTER
from .constants import ALIGNMENT_MISMATCH_CHARACTER
from .constants import SECONDS_BETWEEN_UCSC_REQUESTS
from .crispr_target import check_base_match
from .crispr_target import CrisprAlignment
from .crispr_target import CrisprTarget
from .genomic_sequence import GenomicSequence

__all__ = [
    "GenomicSequence",
    "genomic_sequence",
    "SECONDS_BETWEEN_UCSC_REQUESTS",
    "CrisprTarget",
    "CrisprAlignment",
    "check_base_match",
    "ALIGNMENT_MATCH_CHARACTER",
    "ALIGNMENT_MISMATCH_CHARACTER",
    "ALIGNMENT_GAP_CHARACTER",
]
