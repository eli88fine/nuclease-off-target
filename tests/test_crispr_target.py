# -*- coding: utf-8 -*-
from types import SimpleNamespace

from nuclease_off_target import ALIGNMENT_GAP_CHARACTER
from nuclease_off_target import CAS_VARIETIES
from nuclease_off_target import check_base_match
from nuclease_off_target import create_space_in_alignment_between_guide_and_pam
from nuclease_off_target import crispr_target
from nuclease_off_target import CrisprAlignment
from nuclease_off_target import CrisprTarget
from nuclease_off_target import extract_cigar_str_from_result
from nuclease_off_target import find_all_possible_alignments
from nuclease_off_target import GenomicSequence
from nuclease_off_target import sa_cas_off_target_score
from nuclease_off_target import SaCasTarget
from nuclease_off_target import SEPARATION_BETWEEN_GUIDE_AND_PAM
from nuclease_off_target import SpCasTarget
from nuclease_off_target import VERTICAL_ALIGNMENT_DNA_BULGE_CHARACTER
from nuclease_off_target import VERTICAL_ALIGNMENT_MATCH_CHARACTER
from nuclease_off_target import VERTICAL_ALIGNMENT_MISMATCH_CHARACTER
from nuclease_off_target import VERTICAL_ALIGNMENT_RNA_BULGE_CHARACTER
import pytest
from pytest import approx


def test_CrisprTarget_init_converts_str_sequence_to_BioSeq():
    ct = CrisprTarget("GATTCCGTAGACAGACTAGG", "NGG", -3)
    assert ct.sequence == "GATTCCGTAGACAGACTAGGNGG"


def test_SaCasTarget_init_sets_cutsite_and_pam_and_guide():
    expected_guide = "GCAGAACTACACACCAGGGCC"
    ct = SaCasTarget(expected_guide)
    assert ct.guide_target == expected_guide
    assert (
        ct.cut_site_relative_to_pam == CAS_VARIETIES["Sa"]["cut_site_relative_to_pam"]
    )
    assert ct.pam == CAS_VARIETIES["Sa"]["PAM"]


def test_SpCasTarget_init_sets_cutsite_and_pam_and_guide():
    expected_guide = "GTTGCCCCACAGGGCAGTAA"
    ct = SpCasTarget(expected_guide)
    assert ct.guide_target == expected_guide
    assert (
        ct.cut_site_relative_to_pam == CAS_VARIETIES["Sp"]["cut_site_relative_to_pam"]
    )
    assert ct.pam == CAS_VARIETIES["Sp"]["PAM"]


def test_CrisprAlignment_perform_alignment_on_sequence_that_failed():
    # Eli (10/13/20): root cause of failure was identified as a poor regular expression for finding the right edge of the CIGAR string. The regex needed to be updated to contain a `\D` character to ensure all numbers were captured in group 2
    ct = SaCasTarget("GCAGAACTACACACCAGGGCC")
    gs = GenomicSequence(
        "hg19",
        "chr21",
        43290457,
        True,
        "CACATCTCAAGTGCTCAGAAGCCAGCCCTGGTCTGGAGGCTCAGGGAGTTTTTGTTGAAGCAGGAGGAGCATCTTTCCCGCTGGGCACACGGCCTCGTGACAATATTAAAAGCACTCAGGGAAAACCCACTAACGGAAAAAGCAAACGTATTTTCAAAAATATGCGCCGCAGAACCATCCTGAAAAATGCTGTGGGTGGGCTTCAGGGGGCTAGCGGGTGGGTGGAGACGGACTCTGACTAAGTCAGTTTGTGACACCCTGCAGATTCCCACTCTGGAGTCACGACACACACAGGCACGGCAGAGGCCCCGAAAAGTCCTTCAACAAAAATCTCTGCTGAGATGTGTTCTGTCCAGCAGGTGCCCGAATTTATTTGGGCAAAGAGCACCTTTTGCACACGCTGATGAACGTGCTCAAGAAGGATGGTCCACACTGCAGACTCCCGGGGTCAGCAAAGGGGCAGGGCCTGGCTTTGCTACAAGATAGAGGACTAGACACCAGGCCCAGGGTGTGAGTCTGCCTCTGCTGTTTACCAGATGAGAGAGACCTCAGGTGAGTTATGGCCCCTCCAAGGGCCTCAGTTTCCTTCCCTTACAGTAAGAGGGGTGGCCTGGTACTCACCAAGTCCCTTTCAGTTGGCACAACTTCCAAATTTCAGTTCATTGTTACTGTCTCATAAAGTAGTGTTAAAACACAGATACAGACATACAGATAAAACACAGATACAGACATGTAGATAAAACACAGATACAGATGGGTAGACAGAGGTTTTTTACCCAGCTTCCAAGGCTGTGATGAAATGGCACAATACAGCAACCACTCCTCACTGCCTGCTGCCCAGCAGTTAAAAGGCATGTGTTCTACAGGCCGGCAGTAACTCCCTTTCTTTCTCATTTCCCATTAAATACAATGCTTTCTTCCAGCCTCGTCCAGGGGTTCTGCATGCCCCAGGAGTGTTTGGGGGACTCTCTTTGGCCCTCCATTCCCCCCACCCCCACCCC",
    )
    ca = CrisprAlignment(ct, gs)
    ca.perform_alignment()
    assert ca.genomic_sequence.is_positive_strand is True


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
            "GCAGAACTACACACCAGGGCC",
            "NNGRRT",
            "AGCAGGAGAACATCACACACCAGGGCCTGGGGGTGGG",
            (
                "GCAGAACTACACACCAGGGCCNNGRRT",
                VERTICAL_ALIGNMENT_MISMATCH_CHARACTER * 2
                + VERTICAL_ALIGNMENT_MATCH_CHARACTER
                + VERTICAL_ALIGNMENT_MISMATCH_CHARACTER * 2
                + VERTICAL_ALIGNMENT_MATCH_CHARACTER
                + VERTICAL_ALIGNMENT_MISMATCH_CHARACTER * 2
                + VERTICAL_ALIGNMENT_MATCH_CHARACTER * 18
                + VERTICAL_ALIGNMENT_MISMATCH_CHARACTER,
                "AGAACATCACACACCAGGGCCTGGGGG",
            ),
            "prefers mismatches to two DNA bulges",
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


@pytest.mark.parametrize(
    ",".join(
        (
            "test_guide",
            "test_pam",
            "test_cut_site_relative_to_pam",
            "test_genome_start_coord",
            "test_genome_is_positive_strand",
            "test_genome_seq_str",
            "expected_cut_site",
            "test_description",
        )
    ),
    [
        (
            "GTTAGGACTATTAGCGTGAT",
            "NGG",
            -3,
            100,
            True,
            "AGCGTTAGGACTATTAGCGTGATAGGCTA",
            119,
            "exact match, positive strand, 3 extra genome letters on each side",
        ),
        (
            "GTTAGGACTATTAGCGTGAT",
            "NGG",
            -4,
            110,
            True,
            "TAAGCGTTAGGACTATTAGCGTGATAGGCTAGG",
            130,
            "exact match, positive strand, 5 extra genome letters on each side",
        ),
        (
            "GTTAGGACTATTAGCGTGAT",
            "NGG",
            -3,
            100,
            False,
            "TAGCGTTAGGACTATTAGCGTGATAGGCTA",
            108,
            "exact match, negative strand, 4 extra genome letters on 5' and 3 on 3'",
        ),
        (
            "GTTAGGACTATTAGCGTGAT",
            "NNGRRT",
            -3,
            100,
            True,
            "GTACCGTTAGGACGTATTAGCGTGATCAGAGTCT",
            122,
            "DNA bulge upstream of PAM, positive strand, 5 extra genome letters on 5'",
        ),
        (
            "AGTTAGACTATTAGCGTGAT",
            "NGG",
            -3,
            100,
            True,
            "TTCCGGAGTTGACTATTAGCGTGATAGGTA",
            121,
            "RNA bulge upstream of PAM, positive strand, 6 extra genome letters on 5'",
        ),
    ],
)
def test_CrisprAlignment__align_to_genomic_site__cut_site(
    test_guide,
    test_pam,
    test_cut_site_relative_to_pam,
    test_genome_start_coord,
    test_genome_is_positive_strand,
    test_genome_seq_str,
    expected_cut_site,
    test_description,
):
    ct = CrisprTarget(test_guide, test_pam, test_cut_site_relative_to_pam)
    gs = GenomicSequence(
        "hg19",
        "chr21",
        start_coord=test_genome_start_coord,
        is_positive_strand=test_genome_is_positive_strand,
        sequence=test_genome_seq_str,
    )
    ca = CrisprAlignment(ct, gs)
    ca.perform_alignment()
    assert ca.cut_site_coord == expected_cut_site


@pytest.mark.parametrize(
    ",".join(
        (
            "test_crispr_alignment",
            "test_genome_alignment",
            "expected_score",
            "test_description",
        )
    ),
    [
        (
            "GTTAGGACTATTAGCGTGATNNGRRT",
            "GTTAGGACTATTAGCGTGATAAGAGT",
            0,
            "exact match",
        ),
        (
            "GTTAGGACTATTAGCGTGATNNGRRT",
            "GTTAGGACTATTAGCGTGATAAGAGA",
            2,
            "mismatch in PAM T",
        ),
        (
            "GTTAGGACTATTAGCGTGATNNGRRT",
            "GTTAGGACTATTAGCGTGATAAGACA",
            22,
            "mismatch in PAM T and last R",
        ),
        (
            "GTTAGGACTATTAGCGTGATNNGRRT",
            "GTTAGGACTATTAGCGTGATAAGA" + ALIGNMENT_GAP_CHARACTER + "T",
            20.3,
            "RNA bulge for last R",
        ),
        (
            "GTTAGGACTATTAGCGTGATNNGRRT",
            "GTTAGGACTATTAGCGTGATAAG" + ALIGNMENT_GAP_CHARACTER + "AT",
            20.3,
            "RNA bulge for first R",
        ),
        (
            "GTTAGGACTATTAGCGTGATNNGRRT",
            "GTTAGGACTATTAGCGTGATAA" + ALIGNMENT_GAP_CHARACTER + "AAT",
            40.3,
            "RNA bulge for PAM G",
        ),
        (
            "GTTAGGACTATTAGCGTGATNNGRR" + ALIGNMENT_GAP_CHARACTER + "T",
            "GTTAGGACTATTAGCGTGATAAGAGTT",
            20.3,
            "DNA bulge for last R",
        ),
        (
            "GTTAGGACTATTAGCGTGATNNGR" + ALIGNMENT_GAP_CHARACTER + "RT",
            "GTTAGGACTATTAGCGTGATAAGACGT",
            20.3,
            "DNA bulge for first R",
        ),
        (
            "GTTAGGACTATTAGCGTGATNNG" + ALIGNMENT_GAP_CHARACTER + "RRT",
            "GTTAGGACTATTAGCGTGATAAGTGGT",
            40.3,
            "DNA bulge for PAM G",
        ),
        (
            "GTTAGGACTATTAGCGTGATNN" + ALIGNMENT_GAP_CHARACTER + "GRRT",
            "GTTAGGACTATTAGCGTGATAATGGGT",
            40.3,
            "DNA bulge for last PAM N",
        ),
        (
            "GTTAGGACTATTAGCGTGATN" + ALIGNMENT_GAP_CHARACTER + "NGRRT",
            "GTTAGGACTATTAGCGTGATAATGGGT",
            40.3,
            "DNA bulge for first PAM N",
        ),
        (
            "GTTAGGACTATTAGCGTGATNNGRRT",
            "GTTAGGACTATTAGCGTGATAAGCAA",
            22,
            "mismatch in PAM T and first R",
        ),
        (
            "GTTAGGACTATTAGCGTGATNNGRRT",
            "GTTAGGACTATTAGCGTGATAAGCCT",
            40,
            "mismatch in both PAM Rs",
        ),
        (
            "GTTAGGACTATTAGCGTGATNNGRRT",
            "GTTAGGACTATTAGCGTGATAAAAGT",
            40,
            "mismatch in PAM G",
        ),
        (
            "GTTAGGACTATTAGCGTGATNNGRR" + ALIGNMENT_GAP_CHARACTER + "T",
            "GTTAGGACTATTAGCGTGATAATAGTT",
            60.3,
            "DNA bulge for last R and mismatch in PAM G",
        ),
        (
            "GTTAGGACTATTAGCGTGATNNGRRT",
            "GTTAGGACTATTAGCGTGA" + ALIGNMENT_GAP_CHARACTER + "AAGAAT",
            6.51,
            "RNA bulge just before PAM",
        ),
        (
            "GTTAGGACTATTAGCGTGAT" + ALIGNMENT_GAP_CHARACTER + "NNGRRT",
            "GTTAGGACTATTAGCGTGATCAAGAAT",
            6.7,
            "DNA bulge just before PAM",
        ),
        (
            "GTTAGGACTATTAGCGTGATNNGRRT",
            "GTTAGGACTATTAGCGTG" + ALIGNMENT_GAP_CHARACTER + "TAAGAAT",
            5.51,
            "RNA bulge at position number 2 5' of PAM",
        ),
        (
            "GTTAGGACTATTAGCGTGATNNGRRT",
            "GTTAGGACTATTAGCGTTCTNNGRRT",
            9,
            "Mismatches at positions number 2 and 3 5' of PAM",
        ),
        (
            "GTTAGGACTATTAGCGTGATNNGRRT",
            "GTTAGGACTATTAGCGCGATNNGRRT",
            3,
            "Mismatch at position number 4 5' of PAM",
        ),
        (
            "GTTAGGACTATTAGCGTGATNNGRRT",
            "GCTAGGACTATTAACGTGATNNGRRT",
            1.43,
            "Mismatch at position number 7 and 19 5' of PAM",
        ),
        (
            "ACGTTAGGACTATTAGCGTGATNNGRRT",
            "GGGTTAGGACTATTAGCGTGATNNGRRT",
            0.2,
            "Mismatch at positions number 21 and 22 5' of PAM",
        ),
        (
            "GTTAGGACTATTAGCGTGATNNGRRT",
            "G"
            + ALIGNMENT_GAP_CHARACTER
            + "TAGGACTATTAGCGTGA"
            + ALIGNMENT_GAP_CHARACTER
            + "AAGAAT",
            12.15,
            "RNA bulge just before PAM and near 5' end of guide triggers the extra 5 point 2 bulges penalty",
        ),
        (
            "G" + ALIGNMENT_GAP_CHARACTER + "TTAGGACTATTAGCGTGATNNGRRT",
            "GATTAGGACTATTAGCGTGA" + ALIGNMENT_GAP_CHARACTER + "AAGAAT",
            12.33,
            "RNA bulge just before PAM and DNA bulge near 5' end of guide triggers the extra 5 point 2 bulges penalty",
        ),
    ],
)
def test_sa_cas_off_target_score(
    test_crispr_alignment,
    test_genome_alignment,
    expected_score,
    test_description,
):
    actual = sa_cas_off_target_score((test_crispr_alignment, "", test_genome_alignment))
    assert actual == approx(expected_score)


@pytest.mark.parametrize(
    ",".join(
        (
            "test_crispr_seq",
            "test_genome_seq",
            "test_mismatches",
            "test_total_bulges",
            "test_rna_bulges",
            "test_dna_bulges",
            "expected_alignments",
            "test_description",
        )
    ),
    [
        (
            "GTTAGGACTATTAGCGTGATNNGRRT",
            "GTTAGGACTATTAGCGTGATCCGAGTACG",
            2,
            6,
            3,
            3,
            {
                ("GTTAGGACTATTAGCGTGATNNGRRT", "GTTAGGACTATTAGCGTGATCCGAGT"),
                # --- These two are functionally equivalent, and it may be OK to eliminate in a future version ---
                ("GTTAGGACTATTAGCGTGATNNGRRT", "GTTAGGACTATTAGCGTGATCCG-AG"),
                ("GTTAGGACTATTAGCGTGATNNGRRT", "GTTAGGACTATTAGCGTGATCCGA-G"),
                # -------- #
                ("GTTAGGACTATTAGCGTGATNNGRR-T", "GTTAGGACTATTAGCGTGATCCGAGTA"),
                ("GTTAGGACTATTAGCGTGATNNGRRT", "GTTAGGACTATTAGCGTGATCC-GAG"),
            },
            "finds exact match and some bulges when  bulges allowed",
        ),
        (
            "GTTAGGACTATTAGCGTGATNNGRRT",
            "GTTAGCACTGTTAGCGTGATAAGAGTACG",
            0,
            0,
            0,
            0,
            set(),
            "returns empty set if nothing found",
        ),
        (
            "GCAGAACTACACACCAGGGCCNNGRRT",
            "TTAGAAATACACACTCAGGGCCAGGAATGGTA",
            4,
            1,
            0,
            1,
            {
                (
                    "GCAGAACTACACAC" + ALIGNMENT_GAP_CHARACTER + "CAGGGCCNNGRRT",
                    "TTAGAAATACACACTCAGGGCCAGGAAT",
                )
            },
            "finds DNA bulge",
        ),
        (
            "GCAGAACTACACACCAGGGCCNNGRRT",
            "TGTGAACTACCAGCAGGGCCTAGGGTGGAG",
            5,
            1,
            1,
            0,
            {
                (
                    "GCAGAACTACACACCAGGGCCNNGRRT",
                    "TGTGAACTAC" + ALIGNMENT_GAP_CHARACTER + "CAGCAGGGCCTAGGGT",
                )
            },
            "finds RNA bulge",
        ),
        (
            "GCAGAACTACACACCAGGGCCNNGRRT",
            "GCTTGAACTCAAACCAGGGCCTTGAACATTCCC",
            6,
            2,
            1,
            1,
            {
                ("GCAG-AACTACACACCAGGGCCNNGRRT", "GCTTGAACT-CAAACCAGGGCCTTGAAC"),
                # --- These two are functionally equivalent, and it may be OK to eliminate (likely the one with the more PAM-proximal DNA bulge) in a future version ---
                ("GC-AGAACTACACACCAGGGCCNNGRRT", "GCTTGAACT-CAAACCAGGGCCTTGAAC"),
                ("GCA-GAACTACACACCAGGGCCNNGRRT", "GCTTGAACT-CAAACCAGGGCCTTGAAC"),
                # -------- #
                # --- These two are functionally equivalent, and it may be OK to eliminate (likely the one with the more PAM-proximal DNA bulge) in a future version ---
                ("GCA-GAACTACACACCAGGGCCNNGRRT", "GCTTGAAC-TCAAACCAGGGCCTTGAAC"),
                ("GC-AGAACTACACACCAGGGCCNNGRRT", "GCTTGAAC-TCAAACCAGGGCCTTGAAC"),
                # -------- #
                # --- These two are functionally equivalent, and it may be OK to eliminate (likely the one with the more PAM-proximal DNA bulge) in a future version ---
                ("GCA-GAACTACACACCAGGGCCNNGRRT", "GCTTGAACTC-AAACCAGGGCCTTGAAC"),
                ("GC-AGAACTACACACCAGGGCCNNGRRT", "GCTTGAACTC-AAACCAGGGCCTTGAAC"),
                # -------- #
                # --- These two are functionally equivalent, and it may be OK to eliminate (likely the one with the more PAM-proximal DNA bulge) in a future version ---
                ("GC-AGAACTACACACCAGGGCCNNGRRT", "GCTTGAACTCAA-ACCAGGGCCTTGAAC"),
                ("GCA-GAACTACACACCAGGGCCNNGRRT", "GCTTGAACTCAA-ACCAGGGCCTTGAAC"),
                # -------- #
            },
            "multiple possible alignments",
        ),
        (
            "GCAGAACTACACACCAGGGCCNNGRRT",
            "GCTTGAACTCAAACCAGGGCCTTGAACATTCCC",
            6,
            1,
            1,
            1,
            set(),
            "no alignments found when total bulge limit is less than sum of RNA+DNA bulges",
        ),
        (
            "GCAGAACTACACACCAGGGCCNNGRRT",
            "GCTCGAACTCAAACCAGGGCCTCGAACATTCTCCA",
            5,
            2,
            1,
            1,
            {
                # --- These two are functionally equivalent, and it may be OK to eliminate (likely the one with the more PAM-proximal DNA bulge) in a future version ---
                ("GCA-GAACTACACACCAGGGCCNNGRRT", "GCTCGAACT-CAAACCAGGGCCTCGAAC"),
                ("GC-AGAACTACACACCAGGGCCNNGRRT", "GCTCGAACT-CAAACCAGGGCCTCGAAC"),
                # -------- #
            },
            "another test for multiple possible alignments",
        ),
    ],
)
def test_find_all_possible_alignments(
    test_crispr_seq,
    test_genome_seq,
    test_mismatches,
    test_total_bulges,
    test_rna_bulges,
    test_dna_bulges,
    expected_alignments,
    test_description,
):
    actual_alignments = find_all_possible_alignments(
        test_crispr_seq,
        test_genome_seq,
        test_mismatches,
        test_total_bulges,
        test_rna_bulges,
        test_dna_bulges,
    )
    # print(len(actual_alignments))
    # for al1, al2 in actual_alignments:
    #     print(f"('{al1}',\n '{al2}'),\n")
    assert actual_alignments == expected_alignments


def test_CrisprAlignment__find_optimal_alignment__located_on_positive_strand():
    gs = GenomicSequence(
        "hg19",
        "chr11",
        117756114,
        True,
        "ATCCATCCCTTGTGATTTCTGCATCCTAACCGAGGAGTGAGAAAACACGGAATTTCCTAAATGTTTAAAACATTAAAAAATATTAGCTAACCTATGTATTTGTGAGTGCTGATTATTAATGTCCCCAATAAGGTGAAACTACCAAAGGTATCCAGAGTAAGAGAGGGATCAGAATGGGGATGGGACCGGTCAGGGCCGTCTGCCTGAAGAGCCACCCGTGAAGGAGGAGAGGAGGGTTTAGGGTGGGGGCAGGTAAGCTGAGGCACAGAGGTAGGAAGGACTTGCGTGGGACCTGCAGGCTGCTTACTTGGCAGGTTCAGAAAGTTCTACTGGGAAGGAACAGGGGTTTCAAAGTTTCCATCCAAATAAGACGAAGTCCGTTTGTCTTATTTGGTTCGTCCACCAGCAGAGAGGGGAGAGGACTGGCTGGCGCCCAAGTGGGAGGGCCTTTAACACAGCCGTCCTGGGCCCCACTGTGCTGATAAGAAATCTCACCAGGGCCTTGAAGAGGACCAGATCCAGGCACTAAATCAGCGAGGCGGAGGGCTTCCGAACCCGCTGTCCAGATTCTGCTCTCTGCCTTCCCCTCCCCTCCAAGTCCAGCCTCTCTTTGTCATTGGAAGAGGCCTCGGAGACGCCACCCAACAGAGTTTATTTTTAGTTTTCCATTCTCGGAAGAACCACTGGTCATGGAGCAAGCTGTTTAAATCTTCCCAGCTGTCCTGCTGAGGACTGAACAAGAGAAAAGGGATTTTGCTTGTTGTAGGAAGGGTCAAGGTCAGGTAGCAGGAAGAATTCACCCTGTGGTCAGAGTGGGCAAAATGAGCCAGGTTGGTACATTTATTTCTTTGGCAAATAATTCTTGCTCAACAGAGAGATCAACAATGCTGTCTGCAGTTTCCATGCCAAAGAATGAGGGTCTTTTGGAAGGCAGTAGGGCGGTAGCAGAATGAGCGCAAATAATCATTAATTGAATTCCAGGACTTGCCACTTAATTAC",
    )
    ct = SaCasTarget("GCAGAACTACACACCAGGGCC")
    ca = CrisprAlignment(ct, gs)
    ca.find_optimal_alignment(6, 1, 1, 0)
    assert ca.genomic_sequence.is_positive_strand is True
    assert ca.formatted_alignment == (
        "GCAGAACTACACACCAGGGCCNNGRRT",
        VERTICAL_ALIGNMENT_MISMATCH_CHARACTER * 2
        + VERTICAL_ALIGNMENT_MATCH_CHARACTER * 4
        + VERTICAL_ALIGNMENT_MISMATCH_CHARACTER
        + VERTICAL_ALIGNMENT_MATCH_CHARACTER
        + VERTICAL_ALIGNMENT_RNA_BULGE_CHARACTER
        + VERTICAL_ALIGNMENT_MATCH_CHARACTER
        + VERTICAL_ALIGNMENT_MISMATCH_CHARACTER
        + VERTICAL_ALIGNMENT_MATCH_CHARACTER * 15
        + VERTICAL_ALIGNMENT_MISMATCH_CHARACTER,
        "TAAGAAAT" + ALIGNMENT_GAP_CHARACTER + "CTCACCAGGGCCTTGAAG",
    )
    assert ca.cut_site_coord == 117756114 + 498


def test_CrisprAlignment__find_optimal_alignment__located_on_negative_strand__flips_strand():
    gs = GenomicSequence(
        "hg19",
        "chr8",
        126961707,
        True,
        "ACAGAAATGCAAAACTCACTGATTCATATAAAAAAACAGAGTATTGACTGGGCACGGTGGCTCATGCCTGTAATCTCAGCACTTTGGGAGGCCAAGGCGGTCAGATCATTAAGTCAGAAGATCGAGACCATCCTGGCCAACATGGTGAAACACAGTGTCTACTAAAAATACAAAAATTAGCTGAGCATGGTGGCATGTGCTTGTAATCCCAGCTACTTTGGAGGCTGCGGCAGGAGAATCGCTTGAACCTGGGAGGTAGAGATTGCAGTGAGTCGAGATTGCACCACTGTACACCAGCCTGGCAACAGAGTGAGATTCCATCTCAAAAAACAAACAAATAAACAAACAAACAAAGTATCATCATTGCCTCACGTTTATATTTTTACAAGGATAATATTTTTGGTAGATGGGGCAAGTTGAAGTAACCTGTTCAGTATGAGAACTGAAGCATCTTCACCAATATTTGTGTAGGAGACCCATAAGCTACCCTTGGCTCTGGTGTATAATTCTAGGGAAGCTGTGGACTCAATTTCAAGGGAGGCATGGAGGAGAGACCAAGCGGTCATCTGGTGTCTCCAGAGCCCACCTGAGGGGATAAAACATTCAGCCTGTTTCCGGTGAGCCTTTTTAGCCTCAGGCTGTTGCTCTTAGTGCCCCTGGGCTTCTCTAGAGCCTCCCATGAAGGGACAGGACTTAAAGCTCAAGGGATTATAAGGAATTGAGGGGCCAGAGTGGGAGAGGAAAGGTCTGGCAGGAGTGGATGGAAAGGTGGAGGAGATAGAGGGTTGAAAGAGGAAGATCCAAGGATTAAAGGGAGCTGAAGGGAACATGGAAGATGGTAGAAAAGAGAGATTGGGGGATAGGTGATGAGAAGACACAGGTCACAGATTGGGGAAGGATGACTGTCTTTGAAAGAAGCAATGGTCACTGAAAAGTTTTTGAGAGAGAGCATTTTTTTTTTCTGCCTCTTTGTACCAAATGCCCACCACCACTAA",
    )
    ct = SaCasTarget("GCAGAACTACACACCAGGGCC")
    ca = CrisprAlignment(ct, gs)
    ca.find_optimal_alignment(7, 1, 0, 1)
    assert ca.genomic_sequence.is_positive_strand is False
    assert ca.formatted_alignment == (
        "GCAGAACTACACACCAGGGCCNNGRRT",
        VERTICAL_ALIGNMENT_MISMATCH_CHARACTER * 2
        + VERTICAL_ALIGNMENT_MATCH_CHARACTER * 4
        + VERTICAL_ALIGNMENT_MISMATCH_CHARACTER
        + VERTICAL_ALIGNMENT_MATCH_CHARACTER * 2
        + VERTICAL_ALIGNMENT_MISMATCH_CHARACTER
        + VERTICAL_ALIGNMENT_MATCH_CHARACTER * 7
        + VERTICAL_ALIGNMENT_MISMATCH_CHARACTER
        + VERTICAL_ALIGNMENT_MATCH_CHARACTER * 9,
        "CTAGAATTATACACCAGAGCCAAGGGT",
    )
    assert ca.cut_site_coord == 126962200
