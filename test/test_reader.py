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
from os import path
from glob import glob
import numpy as np

import pytest

import bioread
from bioread import reader
logger = reader.logger

DATA_PATH = path.join(path.dirname(path.abspath(__file__)), "data")

PHYSIO_FILES = set(glob(path.join(DATA_PATH, "physio", "physio-*.acq")))
V3_FILES = set(glob(path.join(DATA_PATH, "physio", "physio-3*.acq")))
V4_FILES = set(glob(path.join(DATA_PATH, "physio", "physio-4*.acq")))

UNCOMPRESSED_FILES = set(glob(path.join(DATA_PATH, "physio", "*[0-9].acq")))
COMPRESSED_FILES = set(glob(path.join(DATA_PATH, "physio", "*-c.acq")))


@pytest.fixture(scope="module")
def uncompressed_datafiles():
    return [
        (f, bioread.read(f)) for f in UNCOMPRESSED_FILES
    ]


def test_greatest_common_denominator():
    assert reader.greatest_common_denominator(8, 12) == 4
    assert reader.greatest_common_denominator(0, 8) == 8


def test_least_common_multiple():
    assert reader.least_common_multiple(4) == 4
    assert reader.least_common_multiple(2, 8) == 8
    assert reader.least_common_multiple(8, 2) == 8
    assert reader.least_common_multiple(2, 7) == 14
    assert reader.least_common_multiple(2, 3, 8) == 24


def assert_pattern(dividers, pattern):
    assert np.array_equal(reader.sample_pattern(dividers), pattern)


def test_sample_pattern():
    assert_pattern([1], [0])
    assert_pattern([1, 2], [0, 1, 0])
    assert_pattern([2, 2], [0, 1])
    assert_pattern([1, 4, 2], [0, 1, 2, 0, 0, 2, 0])


def test_byte_pattern():
    sample_pattern = reader.sample_pattern([1, 2])
    sample_lengths = np.array([1, 4])
    bp = reader.byte_pattern(sample_pattern, sample_lengths)
    assert np.array_equal(bp, np.array([0, 1, 1, 1, 1, 0]))


# def test_reads_files_without_error():
#     for filename in PHYSIO_FILES:
#         logger.debug(filename)
#         assert readers.Reader.read_file(filename)


def test_reads_headers_without_error():
    for filename in PHYSIO_FILES:
        logger.debug("Reading {0}".format(filename))
        df = reader.Reader.read_without_data(filename)
        assert df.samples_per_second is not None
