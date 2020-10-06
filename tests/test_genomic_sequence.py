# -*- coding: utf-8 -*-
import datetime
import time

from freezegun import freeze_time
from nuclease_off_target import genomic_sequence
from nuclease_off_target import GenomicSequence
from nuclease_off_target import SECONDS_BETWEEN_UCSC_REQUESTS


def test_GenomicSequence_init_assigns_values_and_has_str_correct():
    gs = GenomicSequence("hg38", "chr2", 10, True, "ACGATGAT")
    assert str(gs) == "<GenomicSequence hg38 chr2:10-17 +>"


def test_GenomicSequence_init_converts_str_sequence_to_BioSeq():
    gs = GenomicSequence("hg19", "chrX", 22, False, "ATTGGCC")
    assert str(gs) == "<GenomicSequence hg19 chrX:22-28 ->"
    assert gs.sequence.reverse_complement() == "GGCCAAT"


def test_GenomicSequence_from_coordinates__gets_sequence_from_ucsc_browser():
    gs = GenomicSequence.from_coordinates("hg19", "chrX", 2000000, 2000215, True)
    assert (
        str(gs.sequence)
        == "GAGGAAGGAAAAGAAAAGGAAGGAAGGAAAGAAAAGGAAAAGAGGGAAGGAAGGAAGAGAGGAAACATAAAAGGAAGGAAGGAAAGAAGGAAAAGAAAAAGAGGGAAGGAAGGAAAAGAGGAAACATAAAAGGAAGGAAGGAAAGAAGGAAAATAAGAGGAAGGAAGGAAGGAAAGAAGGAGTGAAGGAGGGAGGGAGGGAAGGAGGAAGGGAGGG"
    )


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
