# -*- coding: utf-8 -*-
"""Genomic sequences."""
from Bio.Seq import Seq


class GenomicSequence:
    """Basic definition of a genomic sequence."""

    def __init__(
        self,
        genome: str,
        chromosome: str,
        start_coord: int,
        is_positive_strand: bool,
        sequence: str,
    ) -> None:
        self.genome = genome
        self.chromosome = chromosome
        self.start_coord = start_coord
        self.is_positive_strand = is_positive_strand
        self.sequence = Seq(sequence)
        self.end_coord = self.start_coord + len(self.sequence) - 1

    def __str__(self) -> str:
        return f'<{self.__class__.__name__} {self.genome} {self.chromosome}:{self.start_coord}-{self.end_coord} {"+" if self.is_positive_strand else "-"}>'
