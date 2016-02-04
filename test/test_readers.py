# coding: utf8
# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2016 Board of Regents of the University of Wisconsin System
#
# Written Nate Vack <njvack@wisc.edu> with research from John Ollinger
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison
# Project home: http://github.com/njvack/bioread

from __future__ import absolute_import


from bioread import readers


def test_greatest_common_denominator():
    assert readers.greatest_common_denominator(8, 12) == 4
    assert readers.greatest_common_denominator(0, 8) == 8


def test_least_common_multiple():
    assert readers.least_common_multiple(4) == 4
    assert readers.least_common_multiple(2, 8) == 8
    assert readers.least_common_multiple(8, 2) == 8
    assert readers.least_common_multiple(2, 7) == 14


def assert_pattern(dividers, pattern):
    assert list(readers.sample_pattern(dividers)) == list(pattern)


def test_sample_pattern():
    assert_pattern([1], [0])
    assert_pattern([1, 2], [0, 1, 0])
    assert_pattern([2, 2], [0, 1])
    assert_pattern([1, 4, 2], [0, 1, 2, 0, 0, 2, 0])
