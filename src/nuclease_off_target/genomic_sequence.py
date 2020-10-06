# -*- coding: utf-8 -*-
"""Genomic sequences."""
from Bio.Seq import Seq
from bs4 import BeautifulSoup
import requests


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

    @classmethod
    def from_coordinates(
        cls,
        genome: str,
        chromosome: str,
        start_coord: int,
        end_coord: int,
        is_positive_strand: bool,
    ) -> "GenomicSequence":
        """Create a GenomicSequence from the UCSC Browser."""
        url = f"https://genome.ucsc.edu/cgi-bin/hgc?hgsid=909569459_N8as0yXh8yH3IXZZJcwFBa5u6it3&g=htcGetDna2&table=&i=mixed&getDnaPos={chromosome}%3A{start_coord}-{end_coord}&db={genome}&hgSeq.cdsExon=1&hgSeq.padding5=0&hgSeq.padding3=0&hgSeq.casing=upper&boolshad.hgSeq.maskRepeats=0&hgSeq.repMasking=lower&boolshad.hgSeq.revComp=0&submit=get+DNA"
        response = requests.get(url)
        soup = BeautifulSoup(response.text, "html.parser")
        pre_dom_elements = soup.find_all(
            "pre"
        )  # the genomic text is contained within a <pre> tag
        sequence_element = pre_dom_elements[0]
        sequence_element_lines = str(sequence_element).split("\n")
        lines_of_sequence = sequence_element_lines[2:-1]
        sequence = "".join(lines_of_sequence)
        return cls(genome, chromosome, start_coord, is_positive_strand, sequence)

    def __str__(self) -> str:
        return f'<{self.__class__.__name__} {self.genome} {self.chromosome}:{self.start_coord}-{self.end_coord} {"+" if self.is_positive_strand else "-"}>'
