# -*- coding: utf-8 -*-
from nuclease_off_target import ALIGNMENT_GAP_CHARACTER
from nuclease_off_target import check_base_match
from nuclease_off_target import CrisprAlignment
from nuclease_off_target import CrisprTarget
from nuclease_off_target import GenomicSequence
from nuclease_off_target import VERTICAL_ALIGNMENT_GAP_CHARACTER
from nuclease_off_target import VERTICAL_ALIGNMENT_MATCH_CHARACTER
from nuclease_off_target import VERTICAL_ALIGNMENT_MISMATCH_CHARACTER
import pytest


def test_CrisprTarget_init_converts_str_sequence_to_BioSeq():
    ct = CrisprTarget("GATTCCGTAGACAGACTAGG", "NGG", -3)
    assert ct.sequence == "GATTCCGTAGACAGACTAGGNGG"


def test_CrisprAlignment__align_to_genomic_site__when_perfect_alignment_to_same_strand():
    ct = CrisprTarget("GATTCCGTAGACAGACTAGG", "NGG", -3)
    gs = GenomicSequence(
        "hg19",
        "chr1",
        start_coord=10,
        is_positive_strand=True,
        sequence="AGCTGGATTCCGTAGACAGACTAGGTGGACTG",
    )
    ca = CrisprAlignment(ct, gs)
    ca.perform_alignment()
    assert ca.genomic_sequence.is_positive_strand is True
    assert ca.alignment_result.score == 108
    assert ca.alignment_result.cigar.decode == b"5D20=1X2=4D"


def test_CrisprAlignment__align_to_genomic_site__when_perfect_alignment_to_reverse_strand():
    ct = CrisprTarget("GTTAGGACTATTAGCGTGAT", "NGG", -3)
    gs = GenomicSequence(
        "hg19",
        "chr21",
        start_coord=999,
        is_positive_strand=True,
        sequence="CAGATGTCCCATCACGCTAATAGTCCTAACCGGTTTAG",
    )
    ca = CrisprAlignment(ct, gs)
    ca.perform_alignment()
    assert ca.genomic_sequence.is_positive_strand is False
    assert ca.alignment_result.score == 108
    assert ca.alignment_result.cigar.decode == b"8D20=1X2=7D"
    assert ca.formatted_alignment == (
        "GTTAGGACTATTAGCGTGATNGG",
        VERTICAL_ALIGNMENT_MATCH_CHARACTER * 23,
        "GTTAGGACTATTAGCGTGATGGG",
    )


@pytest.mark.parametrize(
    ",".join(
        (
            "test_guide",
            "test_pam",
            "test_genome_seq_str",
            "expected_formatted_alignment",
            "test_description",
        )
    ),
    [
        (
            "GTTAGGACTATTAGCGTGAT",
            "NGG",
            "GGGTAAGGACTATTAGCGTGATGGGGA",
            (
                "GTTAGGACTATTAGCGTGATNGG",
                VERTICAL_ALIGNMENT_MATCH_CHARACTER * 2
                + VERTICAL_ALIGNMENT_MISMATCH_CHARACTER
                + VERTICAL_ALIGNMENT_MATCH_CHARACTER * 20,
                "GTAAGGACTATTAGCGTGATGGG",
            ),
            "one mismatch",
        ),
        (
            "GTTAGGACTATTAGCGTGAT",
            "RGG",
            "TAGGTCGGGACTATTAGCGTGATTGGATT",
            (
                "GTTAGGACTATTAGCGTGATRGG",
                VERTICAL_ALIGNMENT_MATCH_CHARACTER * 2
                + VERTICAL_ALIGNMENT_MISMATCH_CHARACTER * 2
                + VERTICAL_ALIGNMENT_MATCH_CHARACTER * 16
                + VERTICAL_ALIGNMENT_MISMATCH_CHARACTER
                + VERTICAL_ALIGNMENT_MATCH_CHARACTER * 2,
                "GTCGGGACTATTAGCGTGATTGG",
            ),
            "two mismatches in a row and a third in an ambiguous base",
        ),
        (
            "GTTAGGACTATTAGCGTGAT",
            "NGG",
            "GGGTTGGACTATTAGCGTGATGGGGA",
            (
                "GTTAGGACTATTAGCGTGATNGG",
                VERTICAL_ALIGNMENT_MATCH_CHARACTER * 3
                + VERTICAL_ALIGNMENT_GAP_CHARACTER
                + VERTICAL_ALIGNMENT_MATCH_CHARACTER * 19,
                "GTT" + ALIGNMENT_GAP_CHARACTER + "GGACTATTAGCGTGATGGG",
            ),
            "one RNA bulge",
        ),
        (
            "GTTAGGACTATTAGCGTGAT",
            "NGG",
            "ACCGTTAGGACGTATTAGCGTGATCGGCT",
            (
                "GTTAGGAC" + ALIGNMENT_GAP_CHARACTER + "TATTAGCGTGATNGG",
                VERTICAL_ALIGNMENT_MATCH_CHARACTER * 8
                + VERTICAL_ALIGNMENT_GAP_CHARACTER
                + VERTICAL_ALIGNMENT_MATCH_CHARACTER * 15,
                "GTTAGGACGTATTAGCGTGATCGG",
            ),
            "one DNA bulge",
        ),
    ],
)
def test_CrisprAlignment__align_to_genomic_site__formatted_alignment(
    test_guide,
    test_pam,
    test_genome_seq_str,
    expected_formatted_alignment,
    test_description,
):
    ct = CrisprTarget(test_guide, test_pam, -3)
    gs = GenomicSequence(
        "hg19",
        "chr21",
        start_coord=999,
        is_positive_strand=True,
        sequence=test_genome_seq_str,
    )
    ca = CrisprAlignment(ct, gs)
    ca.perform_alignment()
    assert ca.formatted_alignment == expected_formatted_alignment


@pytest.mark.parametrize(
    ",".join(("test_first_base", "test_second_base", "expected", "test_description")),
    [
        ("A", "A", True, "two equal bases match"),
        ("A", "C", False, "two different bases don't match"),
        ("N", "C", True, "N as first parameter matches to a C"),
        ("N", "A", True, "N as first parameter matches to a A"),
        ("N", "G", True, "N as first parameter matches to a G"),
        ("N", "T", True, "N as first parameter matches to a T"),
        ("R", "A", True, "puRine matches to A"),
        ("R", "G", True, "puRine matches to G"),
        ("R", "C", False, "puRine does not match to C"),
        ("R", "T", False, "puRine does not match to T"),
    ],
)
def test_check_base_match(
    test_first_base, test_second_base, expected, test_description
):
    actual = check_base_match(test_first_base, test_second_base)
    assert actual is expected
