# -*- coding: utf-8 -*-
"""Docstring."""
from . import crispr_target
from . import genomic_sequence
from .constants import ALIGNMENT_GAP_CHARACTER
from .constants import CAS_VARIETIES
from .constants import SECONDS_BETWEEN_UCSC_REQUESTS
from .constants import SEPARATION_BETWEEN_GUIDE_AND_PAM
from .constants import VERTICAL_ALIGNMENT_DNA_BULGE_CHARACTER
from .constants import VERTICAL_ALIGNMENT_MATCH_CHARACTER
from .constants import VERTICAL_ALIGNMENT_MISMATCH_CHARACTER
from .constants import VERTICAL_ALIGNMENT_RNA_BULGE_CHARACTER
from .crispr_target import check_base_match
from .crispr_target import create_space_in_alignment_between_guide_and_pam
from .crispr_target import CrisprAlignment
from .crispr_target import CrisprTarget
from .crispr_target import extract_cigar_str_from_result
from .crispr_target import sa_cas_off_target_score
from .crispr_target import SaCasTarget
from .genomic_sequence import GenomicSequence

__all__ = [
    "GenomicSequence",
    "genomic_sequence",
    "crispr_target",
    "SECONDS_BETWEEN_UCSC_REQUESTS",
    "CrisprTarget",
    "SaCasTarget",
    "CrisprAlignment",
    "sa_cas_off_target_score",
    "check_base_match",
    "extract_cigar_str_from_result",
    "create_space_in_alignment_between_guide_and_pam",
    "VERTICAL_ALIGNMENT_MATCH_CHARACTER",
    "VERTICAL_ALIGNMENT_MISMATCH_CHARACTER",
    "VERTICAL_ALIGNMENT_DNA_BULGE_CHARACTER",
    "VERTICAL_ALIGNMENT_RNA_BULGE_CHARACTER",
    "ALIGNMENT_GAP_CHARACTER",
    "SEPARATION_BETWEEN_GUIDE_AND_PAM",
    "CAS_VARIETIES",
]
