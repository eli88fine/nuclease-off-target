# -*- coding: utf-8 -*-
import datetime
import os
import time

from freezegun import freeze_time
from nuclease_off_target import create_dict_by_chromosome_from_genes
from nuclease_off_target import ExonCoordinates
from nuclease_off_target import GeneCoordinates
from nuclease_off_target import GeneIsoformCoordinates
from nuclease_off_target import genomic_sequence
from nuclease_off_target import GenomicCoordinates
from nuclease_off_target import GenomicSequence
from nuclease_off_target import IsoformInDifferentChromosomeError
from nuclease_off_target import IsoformInDifferentStrandError
from nuclease_off_target import parse_ucsc_refseq_table_into_gene_coordinates
from nuclease_off_target import SECONDS_BETWEEN_UCSC_REQUESTS
import pytest
from stdlib_utils import get_current_file_abs_directory

PATH_OF_CURRENT_FILE = get_current_file_abs_directory()


def test_GenomicSequence_init_assigns_values_and_has_str_correct():
    gs = GenomicSequence("hg38", "chr2", 10, True, "ACGATGAT")
    assert str(gs) == "<GenomicSequence hg38 chr2:10-17 +>"


def test_GenomicSequence_init_converts_str_sequence_to_BioSeq():
    gs = GenomicSequence("hg19", "chrX", 22, False, "ATTGGCC")
    assert str(gs) == "<GenomicSequence hg19 chrX:22-28 ->"
    assert gs.sequence.reverse_complement() == "GGCCAAT"


@pytest.mark.slow
@pytest.mark.timeout(15)
def test_GenomicSequence_from_coordinates__gets_sequence_from_ucsc_browser():
    gs = GenomicSequence.from_coordinates("hg19", "chrX", 2000000, 2000215, True)
    assert (
        str(gs.sequence)
        == "GAGGAAGGAAAAGAAAAGGAAGGAAGGAAAGAAAAGGAAAAGAGGGAAGGAAGGAAGAGAGGAAACATAAAAGGAAGGAAGGAAAGAAGGAAAAGAAAAAGAGGGAAGGAAGGAAAAGAGGAAACATAAAAGGAAGGAAGGAAAGAAGGAAAATAAGAGGAAGGAAGGAAGGAAAGAAGGAGTGAAGGAGGGAGGGAGGGAAGGAGGAAGGGAGGG"
    )


@pytest.mark.slow
@pytest.mark.timeout(15)
def test_GenomicSequence_from_coordinates__gets_sequence_from_ucsc_browser_from_negative_strand():
    gs = GenomicSequence.from_coordinates("hg38", "chr2", 500000, 500010, False)
    assert str(gs.sequence) == "AACCCGTTGGT"


@freeze_time("2020-10-06 14:03:00")
def test_GenomicSequence_from_coordinates__waits_to_ping_browser(mocker):
    mocker.patch.object(
        genomic_sequence,
        "get_time_of_last_request_to_ucsc_browser",
        autospec=True,
        return_value=datetime.datetime(year=2020, month=10, day=6, hour=14, minute=3),
    )
    mocked_request_ucsc = mocker.patch.object(
        genomic_sequence,
        "request_sequence_from_ucsc",
        autospec=True,
        return_value="ACGT",
    )
    mocked_sleep = mocker.patch.object(time, "sleep", autospec=True)
    mocker.patch.object(
        genomic_sequence, "set_time_of_last_request_to_ucsc_browser", autospec=True,
    )
    gs = GenomicSequence.from_coordinates("hg19", "chrX", 2000000, 2000215, True)
    mocked_request_ucsc.assert_called_once()
    mocked_sleep.assert_called_once_with(SECONDS_BETWEEN_UCSC_REQUESTS)
    assert str(gs.sequence) == "ACGT"


@freeze_time(f"2020-10-06 14:03:{SECONDS_BETWEEN_UCSC_REQUESTS}")
def test_GenomicSequence_from_coordinates__does_not_wait_to_ping_ucsc_if_not_needed(
    mocker,
):
    mocker.patch.object(
        genomic_sequence,
        "get_time_of_last_request_to_ucsc_browser",
        autospec=True,
        return_value=datetime.datetime(year=2020, month=10, day=6, hour=14, minute=3),
    )
    mocked_request_ucsc = mocker.patch.object(
        genomic_sequence,
        "request_sequence_from_ucsc",
        autospec=True,
        return_value="ACGT",
    )
    spied_sleep = mocker.spy(time, "sleep")
    mocker.patch.object(
        genomic_sequence, "set_time_of_last_request_to_ucsc_browser", autospec=True,
    )

    GenomicSequence.from_coordinates("hg19", "chrX", 2000000, 2000215, True)
    mocked_request_ucsc.assert_called_once()
    assert spied_sleep.call_count == 0


@freeze_time(f"2020-10-06 14:03:{SECONDS_BETWEEN_UCSC_REQUESTS}")
def test_GenomicSequence_from_coordinates__sets_the_time_of_the_last_ping_to_ucsc(
    mocker,
):
    expected_time = datetime.datetime(
        year=2020,
        month=10,
        day=6,
        hour=14,
        minute=3,
        second=SECONDS_BETWEEN_UCSC_REQUESTS,
    )
    actual_original_time = genomic_sequence.get_time_of_last_request_to_ucsc_browser()
    assert actual_original_time != expected_time
    mocker.patch.object(time, "sleep", autospec=True)
    mocker.patch.object(
        genomic_sequence,
        "request_sequence_from_ucsc",
        autospec=True,
        return_value="ACGT",
    )
    GenomicSequence.from_coordinates("hg19", "chrX", 2000000, 2000215, True)
    actual_new_time = genomic_sequence.get_time_of_last_request_to_ucsc_browser()
    assert actual_new_time == expected_time

    # set it back to the value that it was at the beginning of the test
    genomic_sequence.set_time_of_last_request_to_ucsc_browser(actual_original_time)


def test_GenomicSequence_create_reverse_complement__reverses_an_existing_sequence():
    expected_genome = "hg19"
    expected_chr = "chr2"
    expected_start = 9291912
    expected_stop = 9291940
    gs = GenomicSequence(
        expected_genome, expected_chr, 9291912, False, "TCTTAGACTAGGACGCCGTATGAGCGCAT"
    )
    revcomp = gs.create_reverse_complement()
    assert revcomp.genome == expected_genome
    assert revcomp.chromosome == expected_chr
    assert revcomp.start_coord == expected_start
    assert revcomp.end_coord == expected_stop
    assert revcomp.is_positive_strand is True
    assert str(revcomp.sequence) == "ATGCGCTCATACGGCGTCCTAGTCTAAGA"


def test_GenomicSequence_create_three_prime_trim__on_positive_strand():
    expected_genome = "hg19"
    expected_chr = "chr2"
    expected_start = 9291912
    initial_seq = "TCTTAGACTAGGACGCCGTATGAGCGCAT"
    initial_end = expected_start + len(initial_seq) - 1
    expected_is_positive_strand = True
    gs = GenomicSequence(
        expected_genome,
        expected_chr,
        expected_start,
        expected_is_positive_strand,
        initial_seq,
    )

    # confirm pre-condition
    assert gs.end_coord == initial_end

    actual = gs.create_three_prime_trim(10)
    assert actual != gs
    assert actual.genome == expected_genome
    assert actual.chromosome == expected_chr
    assert actual.start_coord == expected_start
    assert actual.end_coord == initial_end - 10
    assert actual.is_positive_strand is expected_is_positive_strand
    assert str(actual.sequence) == "TCTTAGACTAGGACGCCGT"


def test_GenomicSequence_create_three_prime_trim__on_negative_strand():
    expected_genome = "hg38"
    expected_chr = "chr5"
    expected_start = 188232
    initial_seq = "ACGATGAGAGCCTATATGAGTCTATTGCAGCAGCAGT"
    initial_end = expected_start + len(initial_seq) - 1
    expected_is_positive_strand = False
    gs = GenomicSequence(
        expected_genome,
        expected_chr,
        expected_start,
        expected_is_positive_strand,
        initial_seq,
    )

    # confirm pre-condition
    assert gs.end_coord == initial_end

    actual = gs.create_three_prime_trim(5)
    assert actual != gs
    assert actual.genome == expected_genome
    assert actual.chromosome == expected_chr
    assert actual.start_coord == expected_start + 5
    assert actual.end_coord == initial_end
    assert actual.is_positive_strand is expected_is_positive_strand
    assert str(actual.sequence) == "ACGATGAGAGCCTATATGAGTCTATTGCAGCA"


def test_GenomicSequence_create_five_prime_trim__on_positive_strand():
    expected_genome = "mm10"
    expected_chr = "chr4"
    expected_start = 1234567
    initial_seq = "TATGGAGAGCCGAT"
    initial_end = expected_start + len(initial_seq) - 1
    expected_is_positive_strand = True
    gs = GenomicSequence(
        expected_genome,
        expected_chr,
        expected_start,
        expected_is_positive_strand,
        initial_seq,
    )

    # confirm pre-condition
    assert gs.end_coord == initial_end

    actual = gs.create_five_prime_trim(3)
    assert actual != gs
    assert actual.genome == expected_genome
    assert actual.chromosome == expected_chr
    assert actual.start_coord == expected_start + 3
    assert actual.end_coord == initial_end
    assert actual.is_positive_strand is expected_is_positive_strand
    assert str(actual.sequence) == "GGAGAGCCGAT"


def test_GenomicSequence_create_five_prime_trim__on_negative_strand():
    expected_genome = "hg18"
    expected_chr = "chrX"
    expected_start = 765432
    initial_seq = "ATTGAGCGCTAGCAGCATGGTA"
    initial_end = expected_start + len(initial_seq) - 1
    expected_is_positive_strand = False
    gs = GenomicSequence(
        expected_genome,
        expected_chr,
        expected_start,
        expected_is_positive_strand,
        initial_seq,
    )

    # confirm pre-condition
    assert gs.end_coord == initial_end

    actual = gs.create_five_prime_trim(4)
    assert actual != gs
    assert actual.genome == expected_genome
    assert actual.chromosome == expected_chr
    assert actual.start_coord == expected_start
    assert actual.end_coord == initial_end - 4
    assert actual.is_positive_strand is expected_is_positive_strand
    assert str(actual.sequence) == "AGCGCTAGCAGCATGGTA"


def test_GenomicCoordinates__returns_true_when_equal():
    gc1 = GenomicCoordinates("hg19", "chr2", 9000, 10000)
    gc2 = GenomicCoordinates("hg19", "chr2", 9000, 10000)
    assert gc1 == gc2


def test_ExonCoordinates__from_coordinate_info():
    ec = ExonCoordinates.from_coordinate_info("hg38", "chr4", 5000, 60000, False)
    assert ec.coordinates == GenomicCoordinates("hg38", "chr4", 5000, 60000)
    assert ec.is_positive_strand is False


def test_ExonCoordinates__returns_true_when_equal():
    gc1 = GenomicCoordinates("hg19", "chr2", 9000, 10000)
    gc2 = GenomicCoordinates("hg19", "chr2", 9000, 10000)

    ec1 = ExonCoordinates(gc1, True)
    ec2 = ExonCoordinates(gc2, True)
    assert ec1 == ec2


@pytest.fixture(scope="function", name="generic_negative_strand_gene_isoform")
def fixture_generic_negative_strand_gene_isoform():
    gic = GeneIsoformCoordinates(
        [
            ExonCoordinates.from_coordinate_info("hg38", "chr4", 61000, 62000, False),
            ExonCoordinates.from_coordinate_info("hg38", "chr4", 5000, 60000, False),
        ]
    )
    yield gic


@pytest.fixture(
    scope="function",
    name="generic_negative_strand_gene_isoform_extended_5prime_shortened_3prime_variant",
)
def fixture_generic_negative_strand_gene_isoform_extended_5prime_shortened_3prime_variant():
    gic = GeneIsoformCoordinates(
        [
            ExonCoordinates.from_coordinate_info("hg38", "chr4", 61000, 64000, False),
            ExonCoordinates.from_coordinate_info("hg38", "chr4", 5500, 60000, False),
        ]
    )
    yield gic


@pytest.fixture(
    scope="function",
    name="generic_negative_strand_gene_isoform_extended_3prime_shortened_5prime_variant",
)
def fixture_generic_negative_strand_gene_isoform_extended_3prime_shortened_5prime_variant():
    gic = GeneIsoformCoordinates(
        [
            ExonCoordinates.from_coordinate_info("hg38", "chr4", 61000, 61100, False),
            ExonCoordinates.from_coordinate_info("hg38", "chr4", 2000, 60000, False),
        ]
    )
    yield gic


@pytest.fixture(scope="function", name="generic_positive_strand_gene_isoform")
def fixture_generic_positive_strand_gene_isoform():
    gic = GeneIsoformCoordinates(
        [
            ExonCoordinates.from_coordinate_info("mm10", "chr2", 100000, 102250, True),
            ExonCoordinates.from_coordinate_info("mm10", "chr2", 200000, 203355, True),
        ]
    )
    yield gic


def test_GeneIsoformCoordinates__get_all_exon_coordinates(
    generic_negative_strand_gene_isoform,
):
    actual_coords = generic_negative_strand_gene_isoform.get_all_exon_coordinates()
    assert len(actual_coords) == 2
    assert actual_coords[1].coordinates.end_coord == 60000


def test_GeneIsoformCoordinates__get_start_coord__when_negative_strand(
    generic_negative_strand_gene_isoform,
):
    actual_start = generic_negative_strand_gene_isoform.get_start_coord()
    assert actual_start == 5000


def test_GeneIsoformCoordinates__get_end_coord__when_negative_strand(
    generic_negative_strand_gene_isoform,
):
    actual_end = generic_negative_strand_gene_isoform.get_end_coord()
    assert actual_end == 62000


def test_GeneIsoformCoordinates__get_start_coord__when_positive_strand(
    generic_positive_strand_gene_isoform,
):
    actual_start = generic_positive_strand_gene_isoform.get_start_coord()
    assert actual_start == 100000


def test_GeneIsoformCoordinates__get_end_coord__when_positive_strand(
    generic_positive_strand_gene_isoform,
):
    actual_end = generic_positive_strand_gene_isoform.get_end_coord()
    assert actual_end == 203355


def test_GeneCoordinates__standard_info_can_be_accessed_after_loading_one_isoform(
    generic_negative_strand_gene_isoform,
):
    gc = GeneCoordinates("CCR5", generic_negative_strand_gene_isoform)
    assert gc.is_positive_strand is False
    assert gc.genome == "hg38"
    assert gc.chromosome == "chr4"
    assert gc.get_start_coord() == 5000
    assert gc.get_end_coord() == 62000


def test_GeneCoordinates__get_start_coord__when_isoforms_loaded_after_init(
    generic_negative_strand_gene_isoform,
    generic_negative_strand_gene_isoform_extended_5prime_shortened_3prime_variant,
    generic_negative_strand_gene_isoform_extended_3prime_shortened_5prime_variant,
):
    gc = GeneCoordinates("HBB", generic_negative_strand_gene_isoform)
    # confirm pre-condition
    assert gc.get_start_coord() == 5000
    assert gc.get_end_coord() == 62000

    gc.add_isoform(
        generic_negative_strand_gene_isoform_extended_5prime_shortened_3prime_variant
    )
    assert gc.get_start_coord() == 5000
    assert gc.get_end_coord() == 64000

    gc.add_isoform(
        generic_negative_strand_gene_isoform_extended_3prime_shortened_5prime_variant
    )
    assert gc.get_start_coord() == 2000
    assert gc.get_end_coord() == 64000


def test_GeneCoordinates__add_isoform__raises_error_if_isoform_in_different_chromosome(
    generic_negative_strand_gene_isoform,
):
    gc = GeneCoordinates("CCR2", generic_negative_strand_gene_isoform)
    iso2 = GeneIsoformCoordinates(
        [
            ExonCoordinates.from_coordinate_info("hg38", "chr5", 61000, 62000, False),
            ExonCoordinates.from_coordinate_info("hg38", "chr5", 5000, 60000, False),
        ]
    )

    with pytest.raises(IsoformInDifferentChromosomeError, match=r"CCR2.*chr4.*chr5"):
        gc.add_isoform(iso2)


def test_GeneCoordinates__add_isoform__raises_error_if_isoform_in_different_strand(
    generic_negative_strand_gene_isoform,
):
    gc = GeneCoordinates("HBG", generic_negative_strand_gene_isoform)
    iso2 = GeneIsoformCoordinates(
        [ExonCoordinates.from_coordinate_info("hg38", "chr4", 45000, 50000, True),]
    )

    with pytest.raises(IsoformInDifferentStrandError, match=r"HBG.*negative.*positive"):
        gc.add_isoform(iso2)


def test_GeneIsoformCoordinates__from_ucsc_refseq_table_row__single_exon_positive_strand():
    row = [
        585,
        "NR_036268.1",
        "chr19",
        "+",
        71972,
        "72110",
        "72110",
        "72110",
        "1",
        "71972,",
        "72110,",
        0,
        "MIR1302-11",
        "none",
        "none",
        "-1,",
    ]
    gi = GeneIsoformCoordinates.from_ucsc_refseq_table_row("hg19", row)
    assert gi.genome == "hg19"
    assert gi.is_positive_strand is True
    assert gi.get_start_coord() == 71972
    assert gi.get_end_coord() == 72110
    assert gi.chromosome == "chr19"
    assert len(gi.get_all_exon_coordinates()) == 1


def test_GeneIsoformCoordinates__from_ucsc_refseq_table_row__many_exons_positive_strand():
    row = [
        1,
        "NM_001301825.1",
        "chr1",
        "+",
        33547778,
        33586132,
        33547850,
        33585783,
        9,
        "33547778,33549554,33557650,33558882,33560148,33562307,33563607,33583502,33585644,",
        "33547955,33549728,33557823,33559017,33560314,33562470,33563780,33583717,33586132,",
        "0",
        "AZIN2",
        "cmpl",
        "cmpl",
        "0,0,0,2,2,0,1,0,2,",
    ]
    gi = GeneIsoformCoordinates.from_ucsc_refseq_table_row("hg38", row)
    assert gi.genome == "hg38"
    assert gi.is_positive_strand is True
    assert gi.get_start_coord() == 33547778
    assert gi.get_end_coord() == 33586132
    assert gi.chromosome == "chr1"
    all_exons = gi.get_all_exon_coordinates()
    assert len(all_exons) == 9
    assert all_exons[1].coordinates.start_coord == 33549554
    assert all_exons[3].coordinates.end_coord == 33559017


def test_GeneCoordinates__instance_with_same_name_is_found_in_set(
    generic_negative_strand_gene_isoform,
):
    gc1 = GeneCoordinates("HBD", generic_negative_strand_gene_isoform)
    gc2 = GeneCoordinates("HBD", generic_negative_strand_gene_isoform)
    # check pre-condition
    assert gc1 is not gc2
    actual_set = set([gc1])
    assert gc2 in actual_set


def test_GeneCoordinates__instance_with_same_name_is_equal(
    generic_negative_strand_gene_isoform,
):
    gc1 = GeneCoordinates("HBD", generic_negative_strand_gene_isoform)
    gc2 = GeneCoordinates("HBD", generic_negative_strand_gene_isoform)
    # check pre-condition
    assert gc1 is not gc2

    assert gc1 != 5
    assert gc1 == gc2


def test_GeneCoordinates__instance_with_different_name_is_not_equal(
    generic_negative_strand_gene_isoform,
):
    gc1 = GeneCoordinates("HBD", generic_negative_strand_gene_isoform)
    gc2 = GeneCoordinates("HBB", generic_negative_strand_gene_isoform)
    assert gc1 != gc2


def test_GeneIsoformCoordinates__from_ucsc_refseq_table_row__single_exon_negative_strand():
    row = [
        960,
        "NR_031634.1",
        "chr20",
        "-",
        "49231172",
        49231322,
        49231322,
        "49231322",
        1,
        "49231172,",
        "49231322,",
        "0",
        "MIR1302-5",
        "none",
        "none",
        "-1,",
    ]
    gi = GeneIsoformCoordinates.from_ucsc_refseq_table_row("mm10", row)
    assert gi.genome == "mm10"
    assert gi.is_positive_strand is False
    assert gi.get_start_coord() == 49231172
    assert gi.get_end_coord() == 49231322
    assert gi.chromosome == "chr20"
    assert len(gi.get_all_exon_coordinates()) == 1


def test_parse_ucsc_refseq_table_into_gene_coordinates__puts_multiple_isoforms_into_same_gene():
    filepath = os.path.join(PATH_OF_CURRENT_FILE, "partial_ucsc_hg19_refseq.tsv")
    actual_dict = parse_ucsc_refseq_table_into_gene_coordinates("hg19", filepath)
    actual_keys = list(actual_dict.keys())

    assert actual_dict[actual_keys[0]].genome == "hg19"
    assert len(actual_keys) == 11
    znf215 = actual_dict["ZNF215"]
    assert len(znf215.get_isoforms()) == 10


def test_create_dict_by_chromosome_from_genes():
    filepath = os.path.join(PATH_OF_CURRENT_FILE, "partial_ucsc_hg19_refseq.tsv")
    genes_dict = parse_ucsc_refseq_table_into_gene_coordinates("hg19", filepath)

    genes = genes_dict.values()
    dict_by_chr = create_dict_by_chromosome_from_genes(genes)
    actual_keys = set(dict_by_chr.keys())
    assert actual_keys == set(["chr11", "chr1", "chr9"])
    assert len(dict_by_chr["chr1"]) == 1
    assert len(dict_by_chr["chr11"]) == 4
