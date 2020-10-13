# -*- coding: utf-8 -*-
import datetime
import time

from freezegun import freeze_time
from nuclease_off_target import genomic_sequence
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


@pytest.mark.timeout(15)
def test_GenomicSequence_from_coordinates__gets_sequence_from_ucsc_browser():
    gs = GenomicSequence.from_coordinates("hg19", "chrX", 2000000, 2000215, True)
    assert (
        str(gs.sequence)
        == "GAGGAAGGAAAAGAAAAGGAAGGAAGGAAAGAAAAGGAAAAGAGGGAAGGAAGGAAGAGAGGAAACATAAAAGGAAGGAAGGAAAGAAGGAAAAGAAAAAGAGGGAAGGAAGGAAAAGAGGAAACATAAAAGGAAGGAAGGAAAGAAGGAAAATAAGAGGAAGGAAGGAAGGAAAGAAGGAGTGAAGGAGGGAGGGAGGGAAGGAGGAAGGGAGGG"
    )


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
