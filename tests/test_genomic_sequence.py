# -*- coding: utf-8 -*-
from nuclease_off_target import GenomicSequence


def test_GenomicSequence_init_assigns_values_and_has_str_correct():
    gs = GenomicSequence("hg38", "chr2", 10, True, "ACGATGAT")
    assert str(gs) == "<GenomicSequence hg38 chr2:10-17 +>"


def test_GenomicSequence_init_converts_str_sequence_to_BioSeq():
    gs = GenomicSequence("hg19", "chrX", 22, False, "ATTGGCC")
    assert str(gs) == "<GenomicSequence hg19 chrX:22-28 ->"
    assert gs.sequence.reverse_complement() == "GGCCAAT"
