"""Tests for file_io/exporter.py."""

from __future__ import annotations

import numpy as np

from file_io.exporter import _finite_or_none


def test_finite_or_none_writes_nan_as_null():
    """k-method speeds a branch does not reach are NaN; JSON gets null, not the NaN token."""
    assert _finite_or_none(np.array([1.5, np.nan, -2.0])) == [1.5, None, -2.0]
