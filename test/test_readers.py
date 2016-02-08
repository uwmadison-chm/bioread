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

import pytest

import bioread
from bioread import readers
logger = readers.logger

DATA_PATH = path.join(path.dirname(path.abspath(__file__)), "data")

PHYSIO_FILES = set(glob(path.join(DATA_PATH, "physio", "physio-*.acq")))
V3_FILES = set(glob(path.join(DATA_PATH, "physio", "physio-3*.acq")))
V4_FILES = set(glob(path.join(DATA_PATH, "physio", "physio-4*.acq")))

UNCOMPRESSED_FILES = set(glob(path.join(DATA_PATH, "physio", "*[0-9].acq")))
COMPRESSED_FILES = set(glob(path.join(DATA_PATH, "physio", "*-c.acq")))


@pytest.fixture(scope="module")
def uncompressed_datafiles():
    return [
        (f, bioread.read_file(f)) for f in UNCOMPRESSED_FILES
    ]


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


def test_reads_files_without_error():
    # TODO: This fails both on compressed and v3 files because shitty
    # It's possibly journal data that's causing problems? Probably not
    # for v3 files though.
    for filename in V4_FILES.intersection(UNCOMPRESSED_FILES):
        logger.debug(filename)
        assert readers.AcqReader.read_file(filename)


def test_reads_headers_without_error():
    # TODO This fails with compressed v3 files because shitty
    for filename in PHYSIO_FILES:
        logger.debug("Reading {0}".format(filename))
        df = readers.AcqReader.read_without_data(filename)
        assert df.samples_per_second is not None
