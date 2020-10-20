# -*- coding: utf-8 -*-
"""Genomic sequences."""
from collections import defaultdict
import csv
from dataclasses import dataclass
import datetime
import time
from typing import Any
from typing import Dict
from typing import List
from typing import Sequence
from typing import Set
from typing import Union

from Bio.Seq import Seq
from bs4 import BeautifulSoup
from immutable_data_validation import validate_int
from immutable_data_validation import validate_str
import requests

from .constants import SECONDS_BETWEEN_UCSC_REQUESTS
from .exceptions import IsoformInDifferentChromosomeError
from .exceptions import IsoformInDifferentStrandError

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


@dataclass
class GenomicCoordinates:
    """Coordinates representing a sequence stretch in a genome.

    Start coordinate should always be a lower value than the end
    coordinate.
    """

    genome: str
    chromosome: str
    start_coord: int
    end_coord: int


@dataclass
class ExonCoordinates:
    """Coordinates representing an exon."""

    coordinates: GenomicCoordinates
    is_positive_strand: bool

    @classmethod
    def from_coordinate_info(
        cls,
        genome: str,
        chromosome: str,
        start_coord: int,
        end_coord: int,
        is_positive_strand: bool,
    ) -> "ExonCoordinates":
        """Instantiate from provided details of coordinates."""
        coordinates = GenomicCoordinates(genome, chromosome, start_coord, end_coord)
        return cls(coordinates, is_positive_strand)


class GeneIsoformCoordinates:
    """Coordinates for a specific isoform of a gene."""

    def __init__(self, all_exon_coordinates: Sequence[ExonCoordinates]) -> None:
        self._exon_coordinates = all_exon_coordinates
        self.is_positive_strand = all_exon_coordinates[0].is_positive_strand
        self.genome = all_exon_coordinates[0].coordinates.genome
        self.chromosome = all_exon_coordinates[0].coordinates.chromosome

        self._start_coord: int
        self._end_coord: int
        self._calculate_coordinates_from_exons()

    @classmethod
    def from_ucsc_refseq_table_row(
        cls, genome: str, table_row: Sequence[Union[str, int]]
    ) -> "GeneIsoformCoordinates":
        """Create an instance from a table row.

        Intended to be used on data gathered from the USCS Genome Table Browser using the "RefSeq All (ncbiRefSeq)" table under Track "NCBI RefSeq"
        Example: https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=923625121_aiwBounEVv5j3SwEeuFGRaRYYCOu&clade=mammal&org=&db=hg19&hgta_group=genes&hgta_track=refSeqComposite&hgta_table=ncbiRefSeq&hgta_regionType=genome&position=&hgta_outputType=primaryTable&hgta_outFileName=
        """
        all_exons: List[ExonCoordinates] = list()
        num_exons = validate_int(table_row[8])
        exon_starts = validate_str(table_row[9]).split(",")
        exon_ends = validate_str(table_row[10]).split(",")
        is_positive_strand = table_row[3] == "+"
        chromosome = validate_str(table_row[2])
        for exon_idx in range(num_exons):
            start_coord = validate_int(exon_starts[exon_idx])
            end_coord = validate_int(exon_ends[exon_idx])
            all_exons.append(
                ExonCoordinates.from_coordinate_info(
                    genome, chromosome, start_coord, end_coord, is_positive_strand
                )
            )
        return cls(all_exons)

    def _calculate_coordinates_from_exons(self) -> None:
        """Help init calculate internal attributes."""
        self._start_coord = 10 ** 12  # way larger than any chromosome coordinate
        self._end_coord = 0
        for iter_exon_coords in self._exon_coordinates:
            iter_exon_start = iter_exon_coords.coordinates.start_coord
            iter_exon_end = iter_exon_coords.coordinates.end_coord
            if iter_exon_start < self._start_coord:
                self._start_coord = iter_exon_start
            if iter_exon_end > self._end_coord:
                self._end_coord = iter_exon_end

    def get_all_exon_coordinates(self) -> Sequence[ExonCoordinates]:
        return self._exon_coordinates

    def get_start_coord(self) -> int:
        return self._start_coord

    def get_end_coord(self) -> int:
        return self._end_coord


class GeneCoordinates:
    """Coordinates representing a gene with potentially several isoforms."""

    def __init__(self, name: str, an_isoform: GeneIsoformCoordinates) -> None:
        self.name = name
        self.is_positive_strand = an_isoform.is_positive_strand
        self.genome = an_isoform.genome
        self.chromosome = an_isoform.chromosome
        self._isoforms: Set[GeneIsoformCoordinates] = set()
        self._end_coord = 0
        self._start_coord = 10 ** 12
        self.add_isoform(an_isoform)

    def __hash__(self) -> int:
        return hash(self.name)

    def __eq__(self, other: Any) -> bool:
        if not isinstance(other, self.__class__):
            return False
        return self.name == other.name

    def get_isoforms(self) -> Set[GeneIsoformCoordinates]:
        return self._isoforms

    def add_isoform(self, isoform: GeneIsoformCoordinates) -> None:
        """Add an isoform.

        Perform various validations and calculations.
        """
        if isoform.chromosome != self.chromosome:
            raise IsoformInDifferentChromosomeError(
                f"The current isoforms in this gene ({self.name}) are in chromosome {self.chromosome} but the new isoform attempting to be added is in chromosome {isoform.chromosome}."
            )
        if isoform.is_positive_strand != self.is_positive_strand:
            raise IsoformInDifferentStrandError(
                f"The current isoforms in this gene ({self.name}) are on the {'positive' if self.is_positive_strand else 'negative'} strand but the new isoform attempting to be added is on the {'positive' if isoform.is_positive_strand else 'negative'} strand."
            )
        self._isoforms.add(isoform)
        end_coord = isoform.get_end_coord()
        start_coord = isoform.get_start_coord()
        if end_coord > self._end_coord:
            self._end_coord = end_coord
        if start_coord < self._start_coord:
            self._start_coord = start_coord

    def get_start_coord(self) -> int:
        """Get the most inclusive span of all isoforms."""
        return self._start_coord

    def get_end_coord(self) -> int:
        """Get the most inclusive span of all isoforms."""
        return self._end_coord


def parse_ucsc_refseq_table_into_gene_coordinates(  # pylint:disable=invalid-name # Eli (10/19/20): I know this is a long name
    genome: str, filepath: str,
) -> Dict[str, GeneCoordinates]:
    """Parse the table from the UCSC Genome browser.

    Intended to be used on data gathered from the USCS Genome Table Browser using the "RefSeq All (ncbiRefSeq)" table under Track "NCBI RefSeq"
    Example: https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=923625121_aiwBounEVv5j3SwEeuFGRaRYYCOu&clade=mammal&org=&db=hg19&hgta_group=genes&hgta_track=refSeqComposite&hgta_table=ncbiRefSeq&hgta_regionType=genome&position=&hgta_outputType=primaryTable&hgta_outFileName=

    Returns: a dictionary with the name of the gene as the key and the GeneCoordinates as the value
    """
    dict_of_genes: Dict[str, GeneCoordinates] = dict()
    with open(filepath, newline="") as csvfile:
        the_reader = csv.reader(csvfile, delimiter="\t")
        for row_idx, row in enumerate(the_reader):
            if row_idx == 0:
                continue  # skip header row
            iter_isoform = GeneIsoformCoordinates.from_ucsc_refseq_table_row(
                genome, row
            )
            gene_name = validate_str(row[12])
            iter_gene = GeneCoordinates(gene_name, iter_isoform)
            if gene_name not in dict_of_genes:
                dict_of_genes[gene_name] = iter_gene
            else:
                dict_of_genes[gene_name].add_isoform(iter_isoform)
    return dict_of_genes


def create_dict_by_chromosome_from_genes(
    genes: Sequence[GeneCoordinates],
) -> Dict[str, Set[GeneCoordinates]]:
    dict_by_chr: Dict[str, Set[GeneCoordinates]] = defaultdict(set)
    for iter_gene in genes:
        iter_chromosome = iter_gene.chromosome
        dict_by_chr[iter_chromosome].add(iter_gene)
    return dict_by_chr


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

    def create_three_prime_trim(self, num_bases_to_trim: int) -> "GenomicSequence":
        """Create an instance with some 3' bases trimmed from the sequence."""
        cls = self.__class__
        if self.is_positive_strand:
            new_sequence = cls(
                self.genome,
                self.chromosome,
                self.start_coord,
                self.is_positive_strand,
                str(self.sequence)[:-num_bases_to_trim],
            )
            return new_sequence
        new_sequence = cls(
            self.genome,
            self.chromosome,
            self.start_coord + num_bases_to_trim,
            self.is_positive_strand,
            str(self.sequence)[:-num_bases_to_trim],
        )
        return new_sequence

    def create_five_prime_trim(self, num_bases_to_trim: int) -> "GenomicSequence":
        """Create an instance with some 5' bases trimmed from the sequence."""
        cls = self.__class__
        if self.is_positive_strand:
            new_sequence = cls(
                self.genome,
                self.chromosome,
                self.start_coord + num_bases_to_trim,
                self.is_positive_strand,
                str(self.sequence)[num_bases_to_trim:],
            )
            return new_sequence
        new_sequence = cls(
            self.genome,
            self.chromosome,
            self.start_coord,
            self.is_positive_strand,
            str(self.sequence)[num_bases_to_trim:],
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
