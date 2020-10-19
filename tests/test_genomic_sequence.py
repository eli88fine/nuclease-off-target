# -*- coding: utf-8 -*-
import datetime
import time

from freezegun import freeze_time
from nuclease_off_target import ExonCoordinates
from nuclease_off_target import GeneCoordinates
from nuclease_off_target import GeneIsoformCoordinates
from nuclease_off_target import genomic_sequence
from nuclease_off_target import GenomicCoordinates
from nuclease_off_target import GenomicSequence
from nuclease_off_target import SECONDS_BETWEEN_UCSC_REQUESTS
import pytest


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
