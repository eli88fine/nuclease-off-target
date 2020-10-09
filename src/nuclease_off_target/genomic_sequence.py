# -*- coding: utf-8 -*-
"""Genomic sequences."""
import datetime
import time

from Bio.Seq import Seq
from bs4 import BeautifulSoup
import requests

from .constants import SECONDS_BETWEEN_UCSC_REQUESTS

time_of_last_request_to_ucsc_browser = datetime.datetime(
    year=2019, month=1, day=1
)  # initialize to a value long ago


def get_time_of_last_request_to_ucsc_browser() -> datetime.datetime:  # pylint:disable=invalid-name # Eli (10/6/20): I know this is long, but it's a specific concept that needs a long name
    global time_of_last_request_to_ucsc_browser  # pylint:disable=global-statement,invalid-name # Eli (10/6/20): this is a deliberate use to set up a global singleton
    return time_of_last_request_to_ucsc_browser


def set_time_of_last_request_to_ucsc_browser(  # pylint:disable=invalid-name # Eli (10/6/20): I know this is long, but it's a specific concept that needs a long name
    new_time: datetime.datetime,
) -> None:
    global time_of_last_request_to_ucsc_browser  # pylint:disable=global-statement # Eli (10/6/20): this is a deliberate use to set up a global singleton
    time_of_last_request_to_ucsc_browser = new_time


def request_sequence_from_ucsc(url: str) -> str:
    response = requests.get(url)
    soup = BeautifulSoup(response.text, "html.parser")
    pre_dom_elements = soup.find_all(
        "pre"
    )  # the genomic text is contained within a <pre> tag
    sequence_element = pre_dom_elements[0]
    sequence_element_lines = str(sequence_element).split("\n")
    lines_of_sequence = sequence_element_lines[2:-1]
    return "".join(lines_of_sequence)


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

    def create_reverse_complement(self) -> "GenomicSequence":
        cls = self.__class__
        new_sequence = cls(
            self.genome,
            self.chromosome,
            self.start_coord,
            not self.is_positive_strand,
            str(self.sequence.reverse_complement()),
        )
        return new_sequence

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
        seconds_since_last_call = (
            datetime.datetime.utcnow() - get_time_of_last_request_to_ucsc_browser()
        ).total_seconds()
        seconds_to_wait = SECONDS_BETWEEN_UCSC_REQUESTS - seconds_since_last_call
        if seconds_to_wait > 0:
            time.sleep(seconds_to_wait)
        url = f"https://genome.ucsc.edu/cgi-bin/hgc?hgsid=909569459_N8as0yXh8yH3IXZZJcwFBa5u6it3&g=htcGetDna2&table=&i=mixed&getDnaPos={chromosome}%3A{start_coord}-{end_coord}&db={genome}&hgSeq.cdsExon=1&hgSeq.padding5=0&hgSeq.padding3=0&hgSeq.casing=upper&boolshad.hgSeq.maskRepeats=0&hgSeq.repMasking=lower{'' if is_positive_strand else '&hgSeq.revComp=on'}&boolshad.hgSeq.revComp=1&submit=get+DNA"
        sequence = request_sequence_from_ucsc(url)
        set_time_of_last_request_to_ucsc_browser(datetime.datetime.utcnow())
        return cls(genome, chromosome, start_coord, is_positive_strand, sequence)

    def __str__(self) -> str:
        return f'<{self.__class__.__name__} {self.genome} {self.chromosome}:{self.start_coord}-{self.end_coord} {"+" if self.is_positive_strand else "-"}>'
