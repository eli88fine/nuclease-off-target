# -*- coding: utf-8 -*-
"""Genomic sequences."""
from Bio.Seq import Seq
import parasail
from parasail.bindings_v2 import Result

from .genomic_sequence import GenomicSequence


def check_base_match(possibly_ambiguous_base: str, other_base: str) -> bool:
    """Return true if match."""
    if possibly_ambiguous_base == "N":
        return True
    if possibly_ambiguous_base == other_base:
        return True
    if possibly_ambiguous_base == "R":
        if other_base in ("A", "G"):
            return True
    return False


class CrisprTarget:  # pylint:disable=too-few-public-methods
    def __init__(
        self, guide_target: str, pam: str, cut_site_relative_to_pam: int
    ) -> None:
        self.guide_target = guide_target
        self.pam = pam
        self.cut_site_relative_to_pam = cut_site_relative_to_pam
        self.sequence = Seq(guide_target + pam)


class CrisprAlignment:  # pylint:disable=too-few-public-methods
    """Create an alignment of CRISPR to the Genome."""

    def __init__(
        self, crispr_target: CrisprTarget, genomic_sequence: GenomicSequence
    ) -> None:
        self.crispr_target = crispr_target
        self.genomic_sequence = genomic_sequence
        self.alignment_result: Result

    def perform_alignment(self) -> None:
        """Align CRISPR to Genome."""
        crispr_str = str(self.crispr_target.sequence)
        forward_result = parasail.sg_dx_trace(
            crispr_str, str(self.genomic_sequence.sequence), 10, 10, parasail.dnafull
        )
        genomic_revcomp = self.genomic_sequence.create_reverse_complement()
        revcomp_result = parasail.sg_dx_trace(
            crispr_str, str(genomic_revcomp.sequence), 10, 10, parasail.dnafull
        )
        if forward_result.score >= revcomp_result.score:
            self.alignment_result = forward_result
        else:
            self.genomic_sequence = genomic_revcomp
            self.alignment_result = revcomp_result
