# -*- coding: utf-8 -*-
from nuclease_off_target import ALIGNMENT_GAP_CHARACTER
from nuclease_off_target import SECONDS_BETWEEN_UCSC_REQUESTS
from nuclease_off_target import VERTICAL_ALIGNMENT_DNA_BULGE_CHARACTER
from nuclease_off_target import VERTICAL_ALIGNMENT_MATCH_CHARACTER
from nuclease_off_target import VERTICAL_ALIGNMENT_MISMATCH_CHARACTER
from nuclease_off_target import VERTICAL_ALIGNMENT_RNA_BULGE_CHARACTER


def ucsc():
    assert SECONDS_BETWEEN_UCSC_REQUESTS == 10


def alignment_display():
    assert VERTICAL_ALIGNMENT_MATCH_CHARACTER == " "
    assert VERTICAL_ALIGNMENT_MISMATCH_CHARACTER == "X"
    assert VERTICAL_ALIGNMENT_DNA_BULGE_CHARACTER == "-"
    assert VERTICAL_ALIGNMENT_RNA_BULGE_CHARACTER == "+"
    assert ALIGNMENT_GAP_CHARACTER == "-"
