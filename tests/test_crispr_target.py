# -*- coding: utf-8 -*-
from types import SimpleNamespace

from nuclease_off_target import ALIGNMENT_GAP_CHARACTER
from nuclease_off_target import check_base_match
from nuclease_off_target import create_space_in_alignment_between_guide_and_pam
from nuclease_off_target import crispr_target
from nuclease_off_target import CrisprAlignment
from nuclease_off_target import CrisprTarget
from nuclease_off_target import extract_cigar_str_from_result
from nuclease_off_target import GenomicSequence
from nuclease_off_target import SEPARATION_BETWEEN_GUIDE_AND_PAM
from nuclease_off_target import VERTICAL_ALIGNMENT_DNA_BULGE_CHARACTER
from nuclease_off_target import VERTICAL_ALIGNMENT_MATCH_CHARACTER
from nuclease_off_target import VERTICAL_ALIGNMENT_MISMATCH_CHARACTER
from nuclease_off_target import VERTICAL_ALIGNMENT_RNA_BULGE_CHARACTER
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
    actual_cigar = extract_cigar_str_from_result(ca.alignment_result)
    assert actual_cigar == "5D20=1X2=4D"


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
    actual_cigar = extract_cigar_str_from_result(ca.alignment_result)
    assert actual_cigar == "8D20=1X2=7D"
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
                VERTICAL_ALIGNMENT_MATCH_CHARACTER
                + VERTICAL_ALIGNMENT_MISMATCH_CHARACTER
                + VERTICAL_ALIGNMENT_MATCH_CHARACTER
                + VERTICAL_ALIGNMENT_MISMATCH_CHARACTER
                + VERTICAL_ALIGNMENT_MATCH_CHARACTER * 19,
                "GGTTGGACTATTAGCGTGATGGG",
            ),
            "prefers two mismatches to one RNA bulge",
        ),
        (
            "AGTTAGACTATTAGCGTGAT",
            "NGG",
            "GGAGTTGACTATTAGCGTGATAGGTA",
            (
                "AGTTAGACTATTAGCGTGATNGG",
                VERTICAL_ALIGNMENT_MATCH_CHARACTER * 4
                + VERTICAL_ALIGNMENT_RNA_BULGE_CHARACTER
                + VERTICAL_ALIGNMENT_MATCH_CHARACTER * 18,
                "AGTT" + ALIGNMENT_GAP_CHARACTER + "GACTATTAGCGTGATAGG",
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
                + VERTICAL_ALIGNMENT_DNA_BULGE_CHARACTER
                + VERTICAL_ALIGNMENT_MATCH_CHARACTER * 15,
                "GTTAGGACGTATTAGCGTGATCGG",
            ),
            "one DNA bulge",
        ),
        (
            "GCAGAACTACACACCAGGGCC",
            "NNGRRT",
            "TGTGAGTCCTACCACCAGGGCCTTGGGTCCGA",
            (
                "GCAGAACTACACACCAGGGCCNNGRRT",
                VERTICAL_ALIGNMENT_MISMATCH_CHARACTER * 2
                + VERTICAL_ALIGNMENT_MATCH_CHARACTER * 2
                + VERTICAL_ALIGNMENT_MISMATCH_CHARACTER * 2
                + VERTICAL_ALIGNMENT_MATCH_CHARACTER * 4
                + VERTICAL_ALIGNMENT_RNA_BULGE_CHARACTER
                + VERTICAL_ALIGNMENT_MATCH_CHARACTER * 10
                + VERTICAL_ALIGNMENT_MATCH_CHARACTER * 6,
                "TGAGTCCTAC" + ALIGNMENT_GAP_CHARACTER + "CACCAGGGCCTTGGGT",
            ),
            "ensuring that mismatch occurs at first base instead of a gap later",
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


def test_extract_cigar_str_from_result__windows(mocker):
    mocker.patch.object(
        crispr_target, "is_system_windows", autospec=True, return_value=True
    )
    expected = "some data"
    stub_result = SimpleNamespace(cigar=SimpleNamespace(decode=expected))
    actual = extract_cigar_str_from_result(stub_result)
    assert actual == expected


def test_extract_cigar_str_from_result__linux(mocker):
    mocker.patch.object(
        crispr_target, "is_system_windows", autospec=True, return_value=False
    )
    stub_result = SimpleNamespace(
        cigar=SimpleNamespace(decode=b"some encoded utf-8 data")
    )
    actual = extract_cigar_str_from_result(stub_result)
    assert actual == "some encoded utf-8 data"


@pytest.mark.parametrize(
    ",".join(
        (
            "test_guide",
            "test_pam",
            "initial_formatted_alignment",
            "expected_formatted_alignment",
            "test_description",
        )
    ),
    [
        (
            "GTTAGGACTATTAGCGTGAT",
            "NGG",
            (
                "GTTAGGACTATTAGCGTGATNGG",
                VERTICAL_ALIGNMENT_MATCH_CHARACTER * 2
                + VERTICAL_ALIGNMENT_MISMATCH_CHARACTER
                + VERTICAL_ALIGNMENT_MATCH_CHARACTER * 20,
                "GTAAGGACTATTAGCGTGATGGG",
            ),
            (
                "GTTAGGACTATTAGCGTGAT" + SEPARATION_BETWEEN_GUIDE_AND_PAM + "NGG",
                VERTICAL_ALIGNMENT_MATCH_CHARACTER * 2
                + VERTICAL_ALIGNMENT_MISMATCH_CHARACTER
                + VERTICAL_ALIGNMENT_MATCH_CHARACTER * 17
                + SEPARATION_BETWEEN_GUIDE_AND_PAM
                + VERTICAL_ALIGNMENT_MATCH_CHARACTER * 3,
                "GTAAGGACTATTAGCGTGAT" + SEPARATION_BETWEEN_GUIDE_AND_PAM + "GGG",
            ),
            "one mismatch in guide, SpCas PAM",
        ),
        (
            "GTTAGGACTATTAGCGTGAT",
            "NGG",
            (
                "GTTAGGACTATTAGCGTGATNGG",
                VERTICAL_ALIGNMENT_MATCH_CHARACTER * 21
                + VERTICAL_ALIGNMENT_RNA_BULGE_CHARACTER
                + VERTICAL_ALIGNMENT_MATCH_CHARACTER * 1,
                "GTTAGGACTATTAGCGTGATA" + ALIGNMENT_GAP_CHARACTER + "G",
            ),
            (
                "GTTAGGACTATTAGCGTGAT" + SEPARATION_BETWEEN_GUIDE_AND_PAM + "NGG",
                VERTICAL_ALIGNMENT_MATCH_CHARACTER * 20
                + SEPARATION_BETWEEN_GUIDE_AND_PAM
                + VERTICAL_ALIGNMENT_MATCH_CHARACTER
                + VERTICAL_ALIGNMENT_RNA_BULGE_CHARACTER
                + VERTICAL_ALIGNMENT_MATCH_CHARACTER * 1,
                "GTTAGGACTATTAGCGTGAT"
                + SEPARATION_BETWEEN_GUIDE_AND_PAM
                + "A"
                + ALIGNMENT_GAP_CHARACTER
                + "G",
            ),
            "RNA bulge in PAM, SpCas PAM",
        ),
        (
            "GTTAGGACTATTAGCGTGAT",
            "NGG",
            (
                "GTTAGGACTATTAGCGTGATNG-G",
                VERTICAL_ALIGNMENT_MATCH_CHARACTER * 22
                + VERTICAL_ALIGNMENT_DNA_BULGE_CHARACTER
                + VERTICAL_ALIGNMENT_MATCH_CHARACTER * 1,
                "GTTAGGACTATTAGCGTGATAGTG",
            ),
            (
                "GTTAGGACTATTAGCGTGAT" + SEPARATION_BETWEEN_GUIDE_AND_PAM + "NG-G",
                VERTICAL_ALIGNMENT_MATCH_CHARACTER * 20
                + SEPARATION_BETWEEN_GUIDE_AND_PAM
                + VERTICAL_ALIGNMENT_MATCH_CHARACTER * 2
                + VERTICAL_ALIGNMENT_DNA_BULGE_CHARACTER
                + VERTICAL_ALIGNMENT_MATCH_CHARACTER * 1,
                "GTTAGGACTATTAGCGTGAT" + SEPARATION_BETWEEN_GUIDE_AND_PAM + "AGTG",
            ),
            "one DNA bulge in PAM, SpCas PAM",
        ),
        (
            "GCAGAACTACACACCAGGGCC",
            "NNGRRT",
            (
                "GCAGAACTACACACCAGGGCCNNGRRT",
                VERTICAL_ALIGNMENT_MISMATCH_CHARACTER * 2
                + VERTICAL_ALIGNMENT_MATCH_CHARACTER * 2
                + VERTICAL_ALIGNMENT_MISMATCH_CHARACTER * 2
                + VERTICAL_ALIGNMENT_MATCH_CHARACTER * 4
                + VERTICAL_ALIGNMENT_RNA_BULGE_CHARACTER
                + VERTICAL_ALIGNMENT_MATCH_CHARACTER * 10
                + VERTICAL_ALIGNMENT_MATCH_CHARACTER * 6,
                "TGAGTCCTAC" + ALIGNMENT_GAP_CHARACTER + "CACCAGGGCCTTGGGT",
            ),
            (
                "GCAGAACTACACACCAGGGCC" + SEPARATION_BETWEEN_GUIDE_AND_PAM + "NNGRRT",
                VERTICAL_ALIGNMENT_MISMATCH_CHARACTER * 2
                + VERTICAL_ALIGNMENT_MATCH_CHARACTER * 2
                + VERTICAL_ALIGNMENT_MISMATCH_CHARACTER * 2
                + VERTICAL_ALIGNMENT_MATCH_CHARACTER * 4
                + VERTICAL_ALIGNMENT_RNA_BULGE_CHARACTER
                + VERTICAL_ALIGNMENT_MATCH_CHARACTER * 10
                + SEPARATION_BETWEEN_GUIDE_AND_PAM
                + VERTICAL_ALIGNMENT_MATCH_CHARACTER * 6,
                "TGAGTCCTAC"
                + ALIGNMENT_GAP_CHARACTER
                + "CACCAGGGCC"
                + SEPARATION_BETWEEN_GUIDE_AND_PAM
                + "TTGGGT",
            ),
            "SaCas PAM",
        ),
    ],
)
def test_create_space_in_alignment_between_guide_and_pam(
    test_guide,
    test_pam,
    initial_formatted_alignment,
    expected_formatted_alignment,
    test_description,
):
    ct = CrisprTarget(test_guide, test_pam, -3)
    actual = create_space_in_alignment_between_guide_and_pam(
        initial_formatted_alignment, ct
    )
    assert actual == expected_formatted_alignment
