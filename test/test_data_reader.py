# coding: utf8
# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2025 Board of Regents of the University of Wisconsin System
#
# Written Nate Vack <njvack@wisc.edu> with research from John Ollinger
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison

import numpy as np

from bioread import data_reader


# Lower-level function tests.
def test_greatest_common_denominator():
    assert data_reader.greatest_common_denominator(8, 12) == 4
    assert data_reader.greatest_common_denominator(0, 8) == 8


def test_least_common_multiple():
    assert data_reader.least_common_multiple(4) == 4
    assert data_reader.least_common_multiple(2, 8) == 8
    assert data_reader.least_common_multiple(8, 2) == 8
    assert data_reader.least_common_multiple(2, 7) == 14
    assert data_reader.least_common_multiple(2, 3, 8) == 24


def assert_pattern(dividers, pattern):
    assert np.array_equal(data_reader.sample_pattern(dividers), pattern)


def test_sample_pattern():
    assert_pattern([1], [0])
    assert_pattern([1, 2], [0, 1, 0])
    assert_pattern([2, 2], [0, 1])
    assert_pattern([1, 4, 2], [0, 1, 2, 0, 0, 2, 0])
