# coding: utf8
# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2020 Board of Regents of the University of Wisconsin System
#
# Written Nate Vack <njvack@wisc.edu> with research from John Ollinger
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison
# Project home: http://github.com/njvack/bioread

from __future__ import absolute_import
from os import path
import numpy as np
import itertools
try:
    from html.parser import HTMLParser
except:
    from HTMLParser import HTMLParser  # noqa

import pytest

import bioread
from bioread import reader
from bioread.reader import Reader
import logging
logger = reader.logger
logger.setLevel(logging.DEBUG)

DATA_PATH = path.join(path.dirname(path.abspath(__file__)), "data")

DATASETS = ['physio', 'nojournal']
# 3.8.1 is special, the compressed / uncompressed data don't agree in the last
# chunk. I think it's a bug in acqknowledge.
BADEND_VERSIONS = ['3.8.1']

# Since our data were backconverted to 3.8.1 and then up from there, they're
# all subject to the same bug that bit that version, and they're all different
# than 4.4.0
ORIG_VERSION = ['4.4.0']

# If you have more versions, send me updates!
NORMAL_VERSIONS = [
    '3.9.1',
    '4.1.0',
    '4.2.0',
    '4.3.0',
]

ALL_VERSIONS = BADEND_VERSIONS + NORMAL_VERSIONS + ORIG_VERSION

TEXT_JOURNAL_VERSIONS = [
    '3.8.1',
    '3.9.1',
    '4.1.0'
]

HTML_JOURNAL_VERSIONS = [
    '4.2.0',
    '4.3.0',
    '4.4.0',
    '5.0.1'
]

# Yes this is ridiculous but it helps
COMPRESSIONS = ['', '-c']


def data_file_names(datasets, versions, compressions):
    for dset, ver, comp in itertools.product(
            datasets, versions, compressions):
        yield data_file_name(dset, ver, comp)


def data_file_name(dataset, version, compression):
    return path.join(DATA_PATH, dataset, '{0}-{1}{2}.acq'.format(
        dataset, version, compression))


# Test to make sure we can read everything. This should throw an exception
# if things fail.
@pytest.mark.parametrize(
    'pathname',
    data_file_names(
        DATASETS,
        BADEND_VERSIONS + NORMAL_VERSIONS,
        COMPRESSIONS))
def test_reading(pathname):
    assert Reader.read(pathname), 'Error reading {0}'.format(pathname)


#@pytest.fixture(scope='module')
def canonical_files():
    versions = ['3.8.1', '4.1.0', '4.4.0']
    return dict(
        (
            (ver, dataset), Reader.read(
                data_file_name(dataset, ver, '-c')).datafile
        )
        for ver, dataset in itertools.product(versions, DATASETS)
    )


@pytest.mark.parametrize(
    'test_file,canon_data',
    itertools.chain(
        zip(
            data_file_names(['physio'], ALL_VERSIONS, COMPRESSIONS),
            itertools.repeat(canonical_files()['4.4.0', 'physio'])
        ),
        zip(
            data_file_names(['nojournal'], ALL_VERSIONS, COMPRESSIONS),
            itertools.repeat(canonical_files()['4.4.0', 'nojournal'])

        )
    )
)
def test_full_pattern_channels_match(test_file, canon_data):
    test_data = bioread.read(test_file)
    slices = full_pattern_slices(canon_data.channels)

    for cch, dch, s, i in zip(
            canon_data.channels,
            test_data.channels,
            slices,
            range(len(canon_data.channels))):
        result = np.array_equal(cch.raw_data[s], dch.raw_data[s])
        assert result, '{0} channel {1} does not match'.format(test_file, i)


@pytest.mark.parametrize(
    'dataset,version',
    itertools.product(DATASETS, NORMAL_VERSIONS + ORIG_VERSION)
)
def test_compresed_uncompressed_channels_match(dataset, version):
    uchan = bioread.read(data_file_name(dataset, version, '')).channels
    cchan = bioread.read(data_file_name(dataset, version, '-c')).channels
    for cch, uch, i in zip(uchan, cchan, range(len(uchan))):
        result = np.array_equal(cch.raw_data, uch.raw_data),
        assert result, 'Mismatch for {0}, {1} channel {2}'.format(
            dataset, version, i)


@pytest.mark.parametrize(
    'dataset,version',
    itertools.product(DATASETS, ALL_VERSIONS)
)
def test_compressed_uncompressed_markers_match(dataset, version):
    um = bioread.read_headers(
        data_file_name(dataset, version, '')).event_markers
    cm = bioread.read_headers(
        data_file_name(dataset, version, '-c')).event_markers
    assert um == cm


@pytest.mark.parametrize(
    'test_file,canon_data',
    zip(
        data_file_names(['physio'], HTML_JOURNAL_VERSIONS, COMPRESSIONS),
        itertools.repeat(canonical_files()['4.4.0', 'physio'])
    )
)
def test_html_journals_match(test_file, canon_data):
    test_data = bioread.read(test_file)
    canon_parser = DataExtractor()
    test_parser = DataExtractor()

    canon_parser.feed(canon_data.journal)
    test_parser.feed(test_data.journal)
    assert canon_parser.content == test_parser.content


@pytest.mark.parametrize(
    'test_file,canon_data',
    zip(
        data_file_names(['physio'], TEXT_JOURNAL_VERSIONS, COMPRESSIONS),
        itertools.repeat(canonical_files()['3.8.1', 'physio'])
    )
)
def test_text_journals_match(test_file, canon_data):
    test_data = bioread.read(test_file)

    test_journal = normalize_line_endings(test_data.journal)
    canon_journal = normalize_line_endings(canon_data.journal)

    assert test_journal == canon_journal

def test_reading_r35_file():
    filename = path.join(DATA_PATH, "misc", "r35_test.acq")
    test_data = bioread.read(filename)  # This will raise an exception on fail
    assert len(test_data.channels) == 2

def test_read_iso_8859_1():
    filename = path.join(DATA_PATH, "misc", "iso_8859_1.acq")
    test_data = bioread.read(filename)  # This will raise an exception on fail

    assert len(test_data.channels)==4



# This is kind of intense for something used by tests -- but the deal is:
# different versions of acqknowledge are treating the versions upconverted
# from 3.8.1 differently in the last, partially-filled pattern.
# The data should be identical between all files in a dataset for the
# filled-pattern parts. So this function returns a slice that'll get you that
# part of the raw data.
def full_pattern_slices(channels):
    pattern = reader.sample_pattern([c.frequency_divider for c in channels])
    point_counts = np.array([c.point_count for c in channels])
    pattern_uses = np.bincount(pattern)
    full_pattern_counts = point_counts - (point_counts % pattern_uses)
    return [slice(count) for count in full_pattern_counts]


def normalize_line_endings(s):
    return s.replace('\r\n', '\n').replace('\r', '\n')


# A little thing that'll let us strip text (in a horrible manner) from html
class DataExtractor(HTMLParser):
    def handle_data(self, data):
        if not hasattr(self, 'content'):
            self.content = ''
        self.content += data


# Lower-level function tests.
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
