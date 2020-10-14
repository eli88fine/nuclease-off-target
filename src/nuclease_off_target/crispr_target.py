# -*- coding: utf-8 -*-
"""Genomic sequences."""
import re
from typing import List
from typing import Optional
from typing import Set
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

OUTER_CIGAR_DELETIONS_REGEX = re.compile(r"(\d+)D.*\D(\d+)D")
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


def _count_mismatches_excluding_gaps(
    possibly_ambiguous_sequence: str, other_seq: str
) -> int:
    count = 0
    for idx, base1 in enumerate(possibly_ambiguous_sequence):
        if base1 == ALIGNMENT_GAP_CHARACTER:
            continue
        base2 = other_seq[idx]
        if base2 == ALIGNMENT_GAP_CHARACTER:
            continue
        if not check_base_match(base1, base2):
            count += 1
    return count


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
    cut_site_relative_to_pam = int(CAS_VARIETIES["Sa"]["cut_site_relative_to_pam"])  # type: ignore # Eli (10/12/20) - not sure how to tell mypy here that this will definitely be an int

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


SA_CAS_PAM_POSITIONS_FOR_BULGES = (1, 2, 3, 4, 5)


def sa_cas_off_target_score(alignment: Tuple[str, str, str]) -> Union[float, int]:
    """Calculate COSMID off-target score for SaCas alignment."""
    score: Union[float, int] = 0
    rev_crispr = "".join(reversed(alignment[0]))
    rev_genome = "".join(reversed(alignment[2]))
    crispr_base_position = 0
    guide_mismatch_penalties = CAS_VARIETIES["Sa"][
        "mismatch-penalties-starting-from-PAM"
    ]
    if not isinstance(guide_mismatch_penalties, dict):
        raise NotImplementedError(
            "The mismatch penalties should always be a dictionary."
        )
    for index, crispr_char in enumerate(rev_crispr):
        genome_char = rev_genome[index]
        is_dna_bulge = crispr_char == ALIGNMENT_GAP_CHARACTER
        is_rna_bulge = genome_char == ALIGNMENT_GAP_CHARACTER
        is_mismatch = is_dna_bulge or is_rna_bulge
        if not is_mismatch:
            is_mismatch = not check_base_match(crispr_char, genome_char)
        if is_mismatch:
            if crispr_base_position == 0:
                score += 2
            elif crispr_base_position in (1, 2):
                score += 20
            elif crispr_base_position in (
                3,
                4,
                5,
            ):  # treat any DNA bulges in the "N"s of the PAM as a bulge at the G
                score += 40
            else:
                score += guide_mismatch_penalties[crispr_base_position - 6]
        if is_rna_bulge:
            if crispr_base_position in SA_CAS_PAM_POSITIONS_FOR_BULGES:
                score += 0.3
            else:
                score += 0.51
        if is_dna_bulge:
            if crispr_base_position in SA_CAS_PAM_POSITIONS_FOR_BULGES:
                score += 0.3
            else:
                score += 0.7
        if not is_dna_bulge:
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
    gap_open = 25
    result = parasail.sg_dx_trace(seq1, seq2, gap_open, 30, parasail.dnafull)
    return result


def _is_superfluous_bulge_alignment(crispr_seq: str, genome_seq: str) -> bool:
    if ALIGNMENT_GAP_CHARACTER not in crispr_seq:
        return False
    if ALIGNMENT_GAP_CHARACTER not in genome_seq:
        return False

    mismatch_count = _count_mismatches_excluding_gaps(crispr_seq, genome_seq)

    idx_of_gap_in_crispr = crispr_seq.find(ALIGNMENT_GAP_CHARACTER)
    idx_of_gap_in_genome = genome_seq.find(ALIGNMENT_GAP_CHARACTER)
    crispr_without_gap = (
        crispr_seq[:idx_of_gap_in_crispr] + crispr_seq[idx_of_gap_in_crispr + 1 :]
    )
    genome_without_gap = (
        genome_seq[:idx_of_gap_in_genome] + genome_seq[idx_of_gap_in_genome + 1 :]
    )
    new_mismatch_count = _count_mismatches_excluding_gaps(
        crispr_without_gap, genome_without_gap
    )
    if new_mismatch_count == mismatch_count:
        return True  # Superfluous because both the RNA and DNA bulge could be removed and have the same number of mismatches
    return _is_superfluous_bulge_alignment(crispr_without_gap, genome_without_gap)


def find_all_possible_alignments(
    remaining_crispr_seq: str,
    remaining_genome_seq: str,
    remaining_allowed_mismatches: int,
    remaining_allowed_total_bulges: int,
    remaining_allowed_rna_bulges: int,
    remaining_allowed_dna_bulges: int,
    five_prime_aligned_crispr_seq: str = "",
    five_prime_aligned_genome_seq: str = "",
    found_alignments: Optional[Set[Tuple[str, str]]] = None,
) -> Set[Tuple[str, str]]:
    """Recursively find all possible potential CRISPR/Genome alignments.

    five_prime_aligned_crispr_seq: This will be an empty string at the beginning of the process. Otherwise it will be the result of what's been done in earlier steps of the recursion (including any relevant alignment gap characters)
    five_prime_aligned_genome_seq: This will be an empty string at the beginning of the process. Otherwise it will be the result of what's been done in earlier steps of the recursion (including any relevant alignment gap characters)
    remaining_crispr_seq: This is still yet to be aligned (3' direction)
    remaining_genome_seq: This is still yet to be aligned (3' direction). There should always be some extra genome sequence available to align to to allow for DNA bulges.
    remaining_allowed_mismatches: How many mismatches (including bulges) have yet to be used. When this drops below zero the recursion ends.
    remaining_allowed_total_bulges: How many RNA+DNA bulges have yet to be used. Recursions to create more bulges cease to be spawned when this reaches zero.
    remaining_allowed_rna_bulges: How many RNA bulges have yet to be used. Recursions to create more RNA bulges cease to be spawned when this reaches zero.
    remaining_allowed_dna_bulges: How many DNA bulges have yet to be used. Recursions to create more DNA bulges cease to be spawned when this reaches zero.

    Returns: A tuple of the CRISPR alignment string and Genome alignment string.
    """
    if found_alignments is None:
        found_alignments = set()
    # End the recursion if the whole CRISPR sequence has been aligned
    if len(remaining_crispr_seq) == 0:
        # prune superfluous alignments where offseting bulges could just be deleted
        if not _is_superfluous_bulge_alignment(
            five_prime_aligned_crispr_seq, five_prime_aligned_genome_seq
        ):
            found_alignments.add(
                (five_prime_aligned_crispr_seq, five_prime_aligned_genome_seq)
            )
        return found_alignments
    next_crispr_char = remaining_crispr_seq[0]
    next_genome_char = remaining_genome_seq[0]
    # pylint:disable=too-many-nested-blocks # Eli (10/12/20): keeping it all together for now until the logic is fully worked out for function arguments. Then will refactor to clean up
    if (
        len(five_prime_aligned_crispr_seq) > 1
    ):  # the first character can't be a DNA or RNA bulge
        if (
            remaining_allowed_mismatches > 0
        ):  # bulges also count as a type of mismatch, so there must be allowed mismatches remaining
            if remaining_allowed_total_bulges > 0:
                # Spawn a recursion with an additional RNA bulge
                if remaining_allowed_rna_bulges > 0:
                    if (
                        len(remaining_crispr_seq) > 1
                    ):  # The final CRISPR charactor cannot be a bulge
                        if (
                            five_prime_aligned_genome_seq[-1] != ALIGNMENT_GAP_CHARACTER
                        ):  # Don't allow back-to-back DNA then RNA bulges
                            if (
                                five_prime_aligned_crispr_seq[-1]
                                != ALIGNMENT_GAP_CHARACTER
                            ):  # Don't allow RNA bulges longer than 1 character
                                find_all_possible_alignments(
                                    remaining_crispr_seq[1:],
                                    remaining_genome_seq,
                                    remaining_allowed_mismatches - 1,
                                    remaining_allowed_total_bulges - 1,
                                    remaining_allowed_rna_bulges - 1,
                                    remaining_allowed_dna_bulges,
                                    five_prime_aligned_crispr_seq=five_prime_aligned_crispr_seq
                                    + next_crispr_char,
                                    five_prime_aligned_genome_seq=five_prime_aligned_genome_seq
                                    + ALIGNMENT_GAP_CHARACTER,
                                    found_alignments=found_alignments,
                                )

                # Spawn a recursion with an additional DNA bulge
                if remaining_allowed_dna_bulges > 0:
                    if (
                        five_prime_aligned_crispr_seq[-1] != ALIGNMENT_GAP_CHARACTER
                    ):  # Don't allow back-to-back RNA then DNA bulges
                        if (
                            five_prime_aligned_genome_seq[-1] != ALIGNMENT_GAP_CHARACTER
                        ):  # Don't allow DNA bulges longer than 1 character
                            if (
                                next_crispr_char != "N"
                            ):  # it makes no sense to bulge at an "N"...anything is possible there

                                find_all_possible_alignments(
                                    remaining_crispr_seq,
                                    remaining_genome_seq[1:],
                                    remaining_allowed_mismatches - 1,
                                    remaining_allowed_total_bulges - 1,
                                    remaining_allowed_rna_bulges,
                                    remaining_allowed_dna_bulges - 1,
                                    five_prime_aligned_crispr_seq=five_prime_aligned_crispr_seq
                                    + ALIGNMENT_GAP_CHARACTER,
                                    five_prime_aligned_genome_seq=five_prime_aligned_genome_seq
                                    + next_genome_char,
                                    found_alignments=found_alignments,
                                )

    if not check_base_match(next_crispr_char, next_genome_char):
        remaining_allowed_mismatches -= 1

    # End the recursion if more mismatches have been detected than the limit.
    if remaining_allowed_mismatches < 0:
        return found_alignments
    five_prime_aligned_crispr_seq += next_crispr_char
    five_prime_aligned_genome_seq += next_genome_char
    remaining_crispr_seq = remaining_crispr_seq[1:]
    remaining_genome_seq = remaining_genome_seq[1:]
    return find_all_possible_alignments(
        remaining_crispr_seq,
        remaining_genome_seq,
        remaining_allowed_mismatches,
        remaining_allowed_total_bulges,
        remaining_allowed_rna_bulges,
        remaining_allowed_dna_bulges,
        five_prime_aligned_crispr_seq=five_prime_aligned_crispr_seq,
        five_prime_aligned_genome_seq=five_prime_aligned_genome_seq,
        found_alignments=found_alignments,
    )


def _create_alignment_string(crispr_seq: str, genome_seq: str) -> str:
    """Create the middle alignment string for vertical display."""
    alignment_str = ""
    for char_idx, crispr_char in enumerate(crispr_seq):
        genome_char = genome_seq[char_idx]
        if genome_char == ALIGNMENT_GAP_CHARACTER:
            alignment_str += VERTICAL_ALIGNMENT_RNA_BULGE_CHARACTER
        elif crispr_char == ALIGNMENT_GAP_CHARACTER:
            alignment_str += VERTICAL_ALIGNMENT_DNA_BULGE_CHARACTER
        elif check_base_match(crispr_char, genome_char):
            alignment_str += VERTICAL_ALIGNMENT_MATCH_CHARACTER
        else:
            alignment_str += VERTICAL_ALIGNMENT_MISMATCH_CHARACTER
    return alignment_str


def _get_best_scoring_alignment(
    set_of_alignments: Set[Tuple[str, str]]
) -> Tuple[Tuple[str, str], Union[float, int]]:
    best_scoring_alignment: Optional[Tuple[str, str]] = None
    best_score: Union[float, int] = 9999999
    for crispr_seq, genome_seq in set_of_alignments:
        iter_score = sa_cas_off_target_score((crispr_seq, "", genome_seq))
        if iter_score < best_score:
            best_score = iter_score
            best_scoring_alignment = (crispr_seq, genome_seq)
    if best_scoring_alignment is None:
        raise NotImplementedError(
            "This should never happen...unless no alignments were provided as an input."
        )
    return best_scoring_alignment, best_score


def _find_all_alignments_across_sequence(
    crispr_target_seq: str,
    genome_seq: str,
    allowed_mismatches: int,
    allowed_total_bulges: int,
    allowed_rna_bulges: int,
    allowed_dna_bulges: int,
) -> Set[Tuple[str, str]]:

    max_possible_length_of_genome_alignment = (  # pylint: disable=invalid-name # Eli (10/13/20): I know it's too long
        len(crispr_target_seq) + allowed_dna_bulges
    )

    current_strand_alignments = set()
    for idx in range(len(genome_seq) - max_possible_length_of_genome_alignment):
        iter_genome_seq = genome_seq[
            idx : idx + max_possible_length_of_genome_alignment
        ]
        iter_alignments = find_all_possible_alignments(
            str(crispr_target_seq),
            iter_genome_seq,
            allowed_mismatches,
            allowed_total_bulges,
            allowed_rna_bulges,
            allowed_dna_bulges,
        )
        current_strand_alignments.update(iter_alignments)
    return current_strand_alignments


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

    def find_optimal_alignment(
        self,
        allowed_mismatches: int,
        allowed_total_bulges: int,
        allowed_rna_bulges: int,
        allowed_dna_bulges: int,
    ) -> None:
        """Align CRISPR to genome.

        Searches through both strands to find the optimal alignment (for
        SaCas9). Revcomps the genomic sequnce if the highest scoring
        alignment was on the reverse strand. Calculates cut site
        coordinate and the formatted alignment.
        """
        current_strand_alignments = _find_all_alignments_across_sequence(
            str(self.crispr_target.sequence),
            str(self.genomic_sequence.sequence),
            allowed_mismatches,
            allowed_total_bulges,
            allowed_rna_bulges,
            allowed_dna_bulges,
        )

        best_scoring_current_strand_alignment: Tuple[  # pylint: disable=invalid-name # Eli (10/13/20): I know it's too long
            str, str
        ] = (
            "",
            "",
        )
        best_current_strand_alignment_score: Union[  # pylint: disable=invalid-name # Eli (10/13/20): I know it's too long
            float, int
        ] = 999999
        if current_strand_alignments:
            (
                best_scoring_current_strand_alignment,  # pylint: disable=invalid-name # Eli (10/13/20): I know it's too long
                best_current_strand_alignment_score,  # pylint: disable=invalid-name # Eli (10/13/20): I know it's too long
            ) = _get_best_scoring_alignment(current_strand_alignments)

        opposite_strand_genomic_sequence = (  # pylint: disable=invalid-name # Eli (10/13/20): I know it's too long
            self.genomic_sequence.create_reverse_complement()
        )
        opposite_strand_alignments = _find_all_alignments_across_sequence(
            str(self.crispr_target.sequence),
            str(opposite_strand_genomic_sequence.sequence),
            allowed_mismatches,
            allowed_total_bulges,
            allowed_rna_bulges,
            allowed_dna_bulges,
        )

        best_scoring_opposite_strand_alignment: Tuple[  # pylint: disable=invalid-name # Eli (10/13/20): I know it's too long
            str, str
        ] = (
            "",
            "",
        )
        best_opposite_strand_alignment_score: Union[  # pylint: disable=invalid-name # Eli (10/13/20): I know it's too long
            float, int
        ] = 999999
        if opposite_strand_alignments:
            (
                best_scoring_opposite_strand_alignment,  # pylint: disable=invalid-name # Eli (10/13/20): I know it's too long
                best_opposite_strand_alignment_score,  # pylint: disable=invalid-name # Eli (10/13/20): I know it's too long
            ) = _get_best_scoring_alignment(opposite_strand_alignments)

        best_scoring_alignment = (
            best_scoring_current_strand_alignment
            if best_current_strand_alignment_score
            < best_opposite_strand_alignment_score
            else best_scoring_opposite_strand_alignment
        )
        if best_scoring_alignment == best_scoring_opposite_strand_alignment:
            self.genomic_sequence = opposite_strand_genomic_sequence

        self.formatted_alignment = (
            best_scoring_alignment[0],
            _create_alignment_string(
                best_scoring_alignment[0], best_scoring_alignment[1]
            ),
            best_scoring_alignment[1],
        )

        cut_site_bases_from_three_prime_end = len(  # pylint: disable=invalid-name
            self.crispr_target.pam
        ) + (self.crispr_target.cut_site_relative_to_pam * -1)
        cut_site_index = _find_index_in_alignment_in_crispr_from_three_prime(
            self.formatted_alignment, cut_site_bases_from_three_prime_end
        )
        idx_of_genome_alignment = str(self.genomic_sequence.sequence).find(
            self.formatted_alignment[2].replace(ALIGNMENT_GAP_CHARACTER, "")
        )
        five_prime_genome_seq = (self.formatted_alignment[2][:cut_site_index]).replace(
            ALIGNMENT_GAP_CHARACTER, ""
        )
        trimmed_genomic_seq = self.genomic_sequence.create_three_prime_trim(
            len(self.genomic_sequence.sequence)
            - (idx_of_genome_alignment + len(five_prime_genome_seq))
        )
        if self.genomic_sequence.is_positive_strand:
            self.cut_site_coord = trimmed_genomic_seq.end_coord
        else:
            self.cut_site_coord = (
                trimmed_genomic_seq.start_coord - 1
            )  # adjust so that the cut coordinate is always on the side of the blunt cut towards the coordinate 1 of the chromosome
        # cut_site_index = _find_index_in_alignment_in_crispr_from_three_prime(
        #     self.formatted_alignment, cut_site_bases_from_three_prime_end
        # )
        # # print(five_prime_genome_seq)
        # five_prime_genome_seq = (self.formatted_alignment[2][:cut_site_index]).replace(
        #     ALIGNMENT_GAP_CHARACTER, ""
        # )
        # # print(five_prime_genome_seq)
        # if self.genomic_sequence.is_positive_strand:
        #     self.cut_site_coord = (
        #         self.genomic_sequence.start_coord
        #         + left_count_to_trim
        #         + len(five_prime_genome_seq)
        #         - 1
        #     )  # subtract 1 to get the base 5' of the cut site
        # else:
        #     self.cut_site_coord = (
        #         self.genomic_sequence.end_coord
        #         - left_count_to_trim
        #         - len(five_prime_genome_seq)
        #     )

    def perform_alignment(self) -> None:  # pylint:disable=too-many-locals
        """Align CRISPR to Genome.

        Deprecated
        """
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
        for iter_num_chars, iter_cigar_element_type in cigar_elements:
            iter_num_chars = int(iter_num_chars)
            if iter_cigar_element_type == "=":
                final_crispr_str += temp_crispr_str[:iter_num_chars]
                temp_crispr_str = temp_crispr_str[iter_num_chars:]
                final_genomic_str += temp_genome_str[:iter_num_chars]
                temp_genome_str = temp_genome_str[iter_num_chars:]
            elif iter_cigar_element_type == "X":
                for _ in range(iter_num_chars):
                    crispr_char = temp_crispr_str[0]
                    genome_char = temp_genome_str[0]
                    final_crispr_str += crispr_char
                    temp_crispr_str = temp_crispr_str[1:]
                    final_genomic_str += genome_char
                    temp_genome_str = temp_genome_str[1:]
            elif iter_cigar_element_type == "I":
                if iter_num_chars != 1:
                    raise NotImplementedError("RNA Bulges should only be length of 1")

                crispr_char = temp_crispr_str[0]
                final_crispr_str += crispr_char
                temp_crispr_str = temp_crispr_str[1:]
                final_genomic_str += ALIGNMENT_GAP_CHARACTER
            elif iter_cigar_element_type == "D":
                if iter_num_chars != 1:
                    raise NotImplementedError("DNA Bulges should only be length of 1")
                genome_char = temp_genome_str[0]
                final_genomic_str += genome_char
                temp_genome_str = temp_genome_str[1:]
                final_crispr_str += ALIGNMENT_GAP_CHARACTER

            else:
                raise NotImplementedError(
                    f"Unrecognized cigar element type: {iter_cigar_element_type}"
                )

        # print (cigar_elements)
        alignment_str = _create_alignment_string(final_crispr_str, final_genomic_str)
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
