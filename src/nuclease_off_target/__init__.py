# -*- coding: utf-8 -*-
"""Docstring."""

try:  # adapted from https://packaging.python.org/guides/single-sourcing-package-version/
    from importlib import metadata
except ImportError:  # pragma: no cover
    # Running on pre-3.8 Python; use importlib-metadata package
    import importlib_metadata as metadata  # type: ignore # Eli (9/1/20): for some reason mypy is giving weird errors for this
__version__: str = metadata.version("change_this_to_name_of_package")  # type: ignore # Eli (9/1/20): for some reason mypy is giving weird errors for this

__all__ = ["__version__"]
