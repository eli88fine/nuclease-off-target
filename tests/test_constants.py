# -*- coding: utf-8 -*-
from nuclease_off_target import ALIGNMENT_GAP_CHARACTER
from nuclease_off_target import ALIGNMENT_MATCH_CHARACTER
from nuclease_off_target import ALIGNMENT_MISMATCH_CHARACTER
from nuclease_off_target import SECONDS_BETWEEN_UCSC_REQUESTS


def ucsc():
    assert SECONDS_BETWEEN_UCSC_REQUESTS == 10


def alignment_display():
    assert ALIGNMENT_MATCH_CHARACTER == "|"
    assert ALIGNMENT_MISMATCH_CHARACTER == "X"
    assert ALIGNMENT_GAP_CHARACTER == " "
