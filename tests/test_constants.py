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
    assert SECONDS_BETWEEN_UCSC_REQUESTS == 3


def test_alignment_display():
    assert VERTICAL_ALIGNMENT_MATCH_CHARACTER == " "
    assert VERTICAL_ALIGNMENT_MISMATCH_CHARACTER == "X"
    assert VERTICAL_ALIGNMENT_DNA_BULGE_CHARACTER == "+"
    assert VERTICAL_ALIGNMENT_RNA_BULGE_CHARACTER == "-"
    assert ALIGNMENT_GAP_CHARACTER == "-"
    assert SEPARATION_BETWEEN_GUIDE_AND_PAM == " "


def test_cas_varieties():
    assert CAS_VARIETIES == {
        "Sa": {
            "PAM": "NNGRRT",
            "cut_site_relative_to_pam": -3,
            "mismatch-penalties-starting-from-PAM": {
                0: 6,
                1: 5,
                2: 4,
                3: 3,
                21: 0.1,
                20: 0.1,
                19: 0.12,
                18: 0.13,
                17: 0.15,
                16: 0.17,
                15: 0.19,
                14: 0.21,
                13: 0.23,
                12: 0.27,
                11: 0.35,
                10: 0.5,
                9: 0.7,
                8: 0.8,
                7: 1.1,
                6: 1.3,
                5: 1.9,
                4: 2.3,
            },
        },
        "Sp": {
            "PAM": "NGG",
            "cut_site_relative_to_pam": -3,
            "mismatch-penalties-starting-from-PAM": {
                0: 6,
                1: 5,
                2: 4,
                3: 3,
                21: 0.1,
                20: 0.1,
                19: 0.12,
                18: 0.13,
                17: 0.15,
                16: 0.17,
                15: 0.19,
                14: 0.21,
                13: 0.23,
                12: 0.27,
                11: 0.35,
                10: 0.5,
                9: 0.7,
                8: 0.8,
                7: 1.1,
                6: 1.3,
                5: 1.9,
                4: 2.3,
            },
        },
    }
