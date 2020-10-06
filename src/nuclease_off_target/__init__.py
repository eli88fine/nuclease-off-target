# -*- coding: utf-8 -*-
"""Docstring."""
from . import genomic_sequence
from .constants import SECONDS_BETWEEN_UCSC_REQUESTS
from .genomic_sequence import GenomicSequence

__all__ = ["GenomicSequence", "genomic_sequence", "SECONDS_BETWEEN_UCSC_REQUESTS"]
