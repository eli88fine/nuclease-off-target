# -*- coding: utf-8 -*-
"""Exceptions."""


class IsoformInDifferentChromosomeError(Exception):
    pass


class IsoformInDifferentStrandError(Exception):
    pass


class DnaRequestGenomeMismatchError(Exception):
    def __init__(self, actual_genome: str, expected_genome: str) -> None:
        super().__init__(
            f"The search expected {expected_genome} but found {actual_genome} in the header. Something is wrong with the URL used for the query."
        )


class UrlNotImplementedForGenomeError(NotImplementedError):
    pass
