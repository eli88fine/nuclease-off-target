# -*- coding: utf-8 -*-
from nuclease_off_target import ALIGNMENT_GAP_CHARACTER
from nuclease_off_target import CAS_VARIETIES
from nuclease_off_target import SECONDS_BETWEEN_UCSC_REQUESTS
from nuclease_off_target import SEPARATION_BETWEEN_GUIDE_AND_PAM
from nuclease_off_target import VERTICAL_ALIGNMENT_DNA_BULGE_CHARACTER
from nuclease_off_target import VERTICAL_ALIGNMENT_MATCH_CHARACTER
from nuclease_off_target import VERTICAL_ALIGNMENT_MISMATCH_CHARACTER
from nuclease_off_target import VERTICAL_ALIGNMENT_RNA_BULGE_CHARACTER


def test_ucsc():
    assert SECONDS_BETWEEN_UCSC_REQUESTS == 10


def test_alignment_display():
    assert VERTICAL_ALIGNMENT_MATCH_CHARACTER == " "
    assert VERTICAL_ALIGNMENT_MISMATCH_CHARACTER == "X"
    assert VERTICAL_ALIGNMENT_DNA_BULGE_CHARACTER == "+"
    assert VERTICAL_ALIGNMENT_RNA_BULGE_CHARACTER == "-"
    assert ALIGNMENT_GAP_CHARACTER == "-"
    assert SEPARATION_BETWEEN_GUIDE_AND_PAM == " "


def test_cas_varieties():
    assert CAS_VARIETIES == {"Sa": {"PAM": "NNGRRT", "cut_site_relative_to_pam": -3}}
