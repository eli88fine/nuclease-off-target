# -*- coding: utf-8 -*-
"""Docstring."""
from . import crispr_target
from . import genomic_sequence
from .constants import ALIGNMENT_GAP_CHARACTER
from .constants import SECONDS_BETWEEN_UCSC_REQUESTS
from .constants import VERTICAL_ALIGNMENT_GAP_CHARACTER
from .constants import VERTICAL_ALIGNMENT_MATCH_CHARACTER
from .constants import VERTICAL_ALIGNMENT_MISMATCH_CHARACTER
from .crispr_target import check_base_match
from .crispr_target import CrisprAlignment
from .crispr_target import CrisprTarget
from .crispr_target import extract_cigar_str_from_result
from .genomic_sequence import GenomicSequence

__all__ = [
    "GenomicSequence",
    "genomic_sequence",
    "crispr_target",
    "SECONDS_BETWEEN_UCSC_REQUESTS",
    "CrisprTarget",
    "CrisprAlignment",
    "check_base_match",
    "extract_cigar_str_from_result",
    "VERTICAL_ALIGNMENT_MATCH_CHARACTER",
    "VERTICAL_ALIGNMENT_MISMATCH_CHARACTER",
    "VERTICAL_ALIGNMENT_GAP_CHARACTER",
    "ALIGNMENT_GAP_CHARACTER",
]
