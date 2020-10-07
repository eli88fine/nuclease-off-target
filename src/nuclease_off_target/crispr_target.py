# -*- coding: utf-8 -*-
"""Genomic sequences."""
from Bio.Seq import Seq


class CrisprTarget:  # pylint:disable=too-few-public-methods
    def __init__(
        self, guide_target: str, pam: str, cut_site_relative_to_pam: int
    ) -> None:
        self.guide_target = guide_target
        self.pam = pam
        self.cut_site_relative_to_pam = cut_site_relative_to_pam
        self.sequence = Seq(guide_target + pam)

    # def align_to_genomic_site(self, genomic_site:GenomicSequence):
    #     #result = global_pairwise_align_nucleotide(DNA(str(self.sequence)), DNA(str(genomic_site.sequence)))
    #     #result = parasail.sg_qx(str(self.sequence), str(genomic_site.sequence), 11, 1, parasail.blosum62)
    #     result = parasail.sg_qx('ACTG', 'TTACTG', 11, 1, parasail.blosum62)
    #     print (result.len_query)
