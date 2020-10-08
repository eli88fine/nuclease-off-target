# -*- coding: utf-8 -*-
"""Genomic sequences."""
import re
from typing import Tuple

from Bio.Seq import Seq
import parasail
from parasail.bindings_v2 import Result

from .constants import ALIGNMENT_MATCH_CHARACTER
from .constants import ALIGNMENT_MISMATCH_CHARACTER
from .genomic_sequence import GenomicSequence

OUTER_CIGAR_DELETIONS_REGEX = re.compile(r"(\d+)D.*(\d+)D")
CIGAR_ELEMENT_REGEX = re.compile(r"(\d+)([\=XDI])")


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
        self.formatted_alignment: Tuple[str, str, str]

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
        cigar = self.alignment_result.cigar.decode.decode("utf-8")
        match = OUTER_CIGAR_DELETIONS_REGEX.match(cigar)
        if match is None:
            raise NotImplementedError("There should always be a match to this RegEx.")
        left_count_to_trim = int(match.group(1))
        right_count_to_trim = int(match.group(2))
        trimmed_genomic_seq = str(self.genomic_sequence.sequence)[
            left_count_to_trim:-right_count_to_trim
        ]
        trimmed_cigar = cigar[
            len(str(left_count_to_trim)) + 1 : -(len(str(right_count_to_trim)) + 1)
        ]
        cigar_elements = CIGAR_ELEMENT_REGEX.findall(trimmed_cigar)
        temp_crispr_str = crispr_str
        temp_genome_str = trimmed_genomic_seq
        final_crispr_str = ""
        final_genomic_str = ""
        alignment_str = ""

        for iter_num_chars, iter_cigar_element_type in cigar_elements:
            iter_num_chars = int(iter_num_chars)
            if iter_cigar_element_type == "=":
                alignment_str += iter_num_chars * ALIGNMENT_MATCH_CHARACTER
                final_crispr_str += temp_crispr_str[:iter_num_chars]
                temp_crispr_str = temp_crispr_str[iter_num_chars:]
                final_genomic_str += temp_genome_str[:iter_num_chars]
                temp_genome_str = temp_genome_str[iter_num_chars:]
            elif iter_cigar_element_type == "X":
                for _ in range(iter_num_chars):
                    crispr_char = temp_crispr_str[0]
                    genome_char = temp_genome_str[0]
                    alignment_char = ALIGNMENT_MISMATCH_CHARACTER
                    if check_base_match(crispr_char, genome_char):
                        alignment_char = ALIGNMENT_MATCH_CHARACTER
                    alignment_str += alignment_char
                    final_crispr_str += crispr_char
                    temp_crispr_str = temp_crispr_str[1:]
                    final_genomic_str += genome_char
                    temp_genome_str = temp_genome_str[1:]
            else:
                raise NotImplementedError(
                    f"Unrecognized cigar element type: {iter_cigar_element_type}"
                )

        # print (cigar_elements)
        self.formatted_alignment = (
            final_crispr_str,
            alignment_str,
            final_genomic_str,
        )
