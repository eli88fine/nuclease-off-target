# -*- coding: utf-8 -*-
from nuclease_off_target import CrisprTarget


def test_CrisprTarget_init_converts_str_sequence_to_BioSeq():
    ct = CrisprTarget("GATTCCGTAGACAGACTAGG", "NGG", -3)
    assert ct.sequence == "GATTCCGTAGACAGACTAGGNGG"


# def nottest_CrisprTarget__align_to_genomic_site__when_perfect_alignment():
#     ct = CrisprTarget("GATTCCGTAGACAGACTAGG", "NGG", -3)
#     gs = GenomicSequence(
#         "hg19",
#         "chr1",
#         start_coord=10,
#         is_positive_strand=True,
#         sequence="GATTCCGTAGACAGACTAGGAGG",
#     )
#     alignment_results = ct.align_to_genomic_site(gs)
#     assert alignment_results == 3
