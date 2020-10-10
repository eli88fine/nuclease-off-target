# -*- coding: utf-8 -*-
"""Genomic sequences."""
import re
from typing import List
from typing import Tuple
from typing import Union

from Bio.Seq import Seq
import parasail
from parasail.bindings_v2 import Result
from stdlib_utils import is_system_windows

from .constants import ALIGNMENT_GAP_CHARACTER
from .constants import CAS_VARIETIES
from .constants import SEPARATION_BETWEEN_GUIDE_AND_PAM
from .constants import VERTICAL_ALIGNMENT_DNA_BULGE_CHARACTER
from .constants import VERTICAL_ALIGNMENT_MATCH_CHARACTER
from .constants import VERTICAL_ALIGNMENT_MISMATCH_CHARACTER
from .constants import VERTICAL_ALIGNMENT_RNA_BULGE_CHARACTER
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


class SaCasTarget(CrisprTarget):  # pylint:disable=too-few-public-methods
    pam = str(CAS_VARIETIES["Sa"]["PAM"])
    cut_site_relative_to_pam = int(CAS_VARIETIES["Sa"]["cut_site_relative_to_pam"])

    def __init__(self, guide_target: str) -> None:
        super().__init__(guide_target, self.pam, self.cut_site_relative_to_pam)


def extract_cigar_str_from_result(result: Result) -> str:
    """Extract the CIGAR alignment from the Parasail alignment Result.

    For some reason, in Windows the result is not bytes-encoded but in
    Linux it is. So needed to have a way to handle both.
    """
    if is_system_windows():
        cigar = result.cigar.decode
        if not isinstance(cigar, str):
            raise NotImplementedError("The decoded CIGAR should always be a string.")
        return cigar
    cigar = result.cigar.decode.decode("utf-8")
    if not isinstance(cigar, str):
        raise NotImplementedError("The decoded CIGAR should always be a string.")
    return cigar


def _find_index_in_alignment_in_crispr_from_three_prime(  # pylint: disable=invalid-name
    alignment: Tuple[str, str, str], num_bases: int
) -> int:
    found_bases_count = 0
    crispr_str = alignment[0]
    for start_idx in range(len(crispr_str) - 1, 0, -1):
        if crispr_str[start_idx] != ALIGNMENT_GAP_CHARACTER:
            found_bases_count += 1
        if found_bases_count == num_bases:
            break
    else:
        raise NotImplementedError(
            "The loop should never complete normally---enough letters to create the PAM and any additional sequence should always be found."
        )
    return start_idx


def sa_cas_off_target_score(alignment: Tuple[str, str, str]) -> Union[float, int]:
    """Calculate COSMID off-target score for SaCas alignment."""
    score = 0
    rev_crispr = "".join(reversed(alignment[0]))
    rev_genome = "".join(reversed(alignment[2]))
    crispr_base_position = 0
    for index, crispr_char in enumerate(rev_crispr):
        genome_char = rev_genome[index]
        is_dna_bulge = crispr_char == ALIGNMENT_GAP_CHARACTER
        # is_rna_bulge = genome_char == ALIGNMENT_GAP_CHARACTER
        is_mismatch = is_dna_bulge
        # if not is_dna_bulge:
        is_mismatch = not check_base_match(crispr_char, genome_char)
        if crispr_base_position == 0:
            if is_mismatch:
                score += 2
        elif crispr_base_position in (1, 2):
            if is_mismatch:
                score += 20
        else:
            score += 0
        # if not is_dna_bulge:
        crispr_base_position += 1
    return score


def create_space_in_alignment_between_guide_and_pam(  # pylint:disable=invalid-name # Eli (10/9/20): I know this is too long, but unsure a better way to describe it
    alignment: Tuple[str, str, str], crispr_target: CrisprTarget
) -> Tuple[str, str, str]:
    """Adjust an alignment to create visual space between Guide and PAM."""
    pam_len = len(crispr_target.pam)
    pam_start_idx = _find_index_in_alignment_in_crispr_from_three_prime(
        alignment, pam_len
    )
    new_alignment: List[str] = list()
    for iter_alignment in alignment:
        new_alignment.append(
            iter_alignment[:pam_start_idx]
            + SEPARATION_BETWEEN_GUIDE_AND_PAM
            + iter_alignment[pam_start_idx:]
        )
    return (new_alignment[0], new_alignment[1], new_alignment[2])


def _run_alignment(seq1: str, seq2: str) -> Result:
    gap_open = 15
    result = parasail.sg_dx_trace(seq1, seq2, gap_open, 10, parasail.dnafull)
    return result


class CrisprAlignment:  # pylint:disable=too-few-public-methods
    """Create an alignment of CRISPR to the Genome."""

    def __init__(
        self, crispr_target: CrisprTarget, genomic_sequence: GenomicSequence
    ) -> None:
        self.crispr_target = crispr_target
        self.genomic_sequence = genomic_sequence
        self.alignment_result: Result
        self.formatted_alignment: Tuple[str, str, str]
        self.cut_site_coord: int  # the base 5' (on positive strand...so always closer to start coordinate of chromosome) of the blunt cut

    def perform_alignment(self) -> None:  # pylint:disable=too-many-locals
        """Align CRISPR to Genome."""
        crispr_str = str(self.crispr_target.sequence)
        forward_result = _run_alignment(crispr_str, str(self.genomic_sequence.sequence))
        genomic_revcomp = self.genomic_sequence.create_reverse_complement()
        revcomp_result = _run_alignment(crispr_str, str(genomic_revcomp.sequence))
        if forward_result.score >= revcomp_result.score:
            self.alignment_result = forward_result
        else:
            self.genomic_sequence = genomic_revcomp
            self.alignment_result = revcomp_result
        cigar = extract_cigar_str_from_result(self.alignment_result)
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
                alignment_str += iter_num_chars * VERTICAL_ALIGNMENT_MATCH_CHARACTER
                final_crispr_str += temp_crispr_str[:iter_num_chars]
                temp_crispr_str = temp_crispr_str[iter_num_chars:]
                final_genomic_str += temp_genome_str[:iter_num_chars]
                temp_genome_str = temp_genome_str[iter_num_chars:]
            elif iter_cigar_element_type == "X":
                for _ in range(iter_num_chars):
                    crispr_char = temp_crispr_str[0]
                    genome_char = temp_genome_str[0]
                    alignment_char = VERTICAL_ALIGNMENT_MISMATCH_CHARACTER
                    if check_base_match(crispr_char, genome_char):
                        alignment_char = VERTICAL_ALIGNMENT_MATCH_CHARACTER
                    alignment_str += alignment_char
                    final_crispr_str += crispr_char
                    temp_crispr_str = temp_crispr_str[1:]
                    final_genomic_str += genome_char
                    temp_genome_str = temp_genome_str[1:]
            elif iter_cigar_element_type == "I":
                if iter_num_chars != 1:
                    raise NotImplementedError("RNA Bulges should only be length of 1")

                crispr_char = temp_crispr_str[0]
                alignment_str += VERTICAL_ALIGNMENT_RNA_BULGE_CHARACTER
                final_crispr_str += crispr_char
                temp_crispr_str = temp_crispr_str[1:]
                final_genomic_str += ALIGNMENT_GAP_CHARACTER
            elif iter_cigar_element_type == "D":
                if iter_num_chars != 1:
                    raise NotImplementedError("DNA Bulges should only be length of 1")
                genome_char = temp_genome_str[0]
                alignment_str += VERTICAL_ALIGNMENT_DNA_BULGE_CHARACTER
                final_genomic_str += genome_char
                temp_genome_str = temp_genome_str[1:]
                final_crispr_str += ALIGNMENT_GAP_CHARACTER

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
        # print("\n")
        # for line in self.formatted_alignment:
        #     print(line)

        cut_site_bases_from_three_prime_end = len(  # pylint: disable=invalid-name
            self.crispr_target.pam
        ) + (self.crispr_target.cut_site_relative_to_pam * -1)
        cut_site_index = _find_index_in_alignment_in_crispr_from_three_prime(
            self.formatted_alignment, cut_site_bases_from_three_prime_end
        )
        # print(five_prime_genome_seq)
        five_prime_genome_seq = (self.formatted_alignment[2][:cut_site_index]).replace(
            ALIGNMENT_GAP_CHARACTER, ""
        )
        # print(five_prime_genome_seq)
        if self.genomic_sequence.is_positive_strand:
            self.cut_site_coord = (
                self.genomic_sequence.start_coord
                + left_count_to_trim
                + len(five_prime_genome_seq)
                - 1
            )  # subtract 1 to get the base 5' of the cut site
        else:
            self.cut_site_coord = (
                self.genomic_sequence.end_coord
                - left_count_to_trim
                - len(five_prime_genome_seq)
            )
