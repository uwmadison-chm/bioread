# coding: utf8
# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2025 Board of Regents of the University of Wisconsin System
#
# Written Nate Vack <njvack@wisc.edu> with research from John Ollinger
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison

import pytest

from os import path
import numpy as np
from html.parser import HTMLParser

import bioread
from bioread import reader
from bioread.reader import Reader

import file_groups


def test_read_does_not_raise_exception(any_acq_file):
    assert Reader.read(any_acq_file) is not None


def test_text_journals_match(text_journal_file, journal_text):
    data = bioread.read(text_journal_file)
    assert data is not None
    assert clean_text(data.journal) == clean_text(journal_text)


def test_html_journals_match(html_journal_file, journal_html):
    data = bioread.read(html_journal_file)
    saved_html_extractor = HTMLTextExtractor()
    saved_html_extractor.feed(journal_html)
    acq_html_extractor = HTMLTextExtractor()
    acq_html_extractor.feed(data.journal)

    assert saved_html_extractor.content == acq_html_extractor.content


def test_biopac_object_smoke(any_acq_file):
    """
    Stringifying events and channels exercises a lot of stuff actually
    """
    reader = Reader.read(any_acq_file)
    assert reader.read_errors == []
    data = reader.datafile
    assert len(data.channels) > 0
    assert len(data.event_markers) > 0
    assert str(data)
    assert [str(ch) for ch in data.channels]
    assert [str(m) for m in data.event_markers]


def test_compressed_uncompressed_good_channels_match(good_data_compressed_pair):
    cfile, ufile = good_data_compressed_pair
    cdata = bioread.read(cfile)
    udata = bioread.read(ufile)
    for cch, uch in zip(cdata.channels, udata.channels):
        assert np.array_equal(cch.raw_data, uch.raw_data)


def test_compressed_uncompressed_buggy_channels_match(buggy_data_compressed_pair):
    cfile, ufile = buggy_data_compressed_pair
    cdata = bioread.read(cfile)
    udata = bioread.read(ufile)
    # Because of the bad end of the uncompressed file, only compare the first
    # 90% of the data
    valid_slices = full_pattern_slices(cdata.channels)
    for cch, uch, channel_slice in zip(cdata.channels, udata.channels, valid_slices):
        assert np.array_equal(cch.raw_data[channel_slice], uch.raw_data[channel_slice])


def test_compressed_uncompressed_markers_match(compressed_uncompressed_pair):
    """
    The compressed 5.0.1 physio file seems to have messed up marker headers, so we
    skip it.
    """
    cfile, ufile = compressed_uncompressed_pair
    if 'physio-5.0.1' in str(cfile):
        pytest.xfail("5.0.1 compressed physio file has marker header issues")
    cm = bioread.read_headers(cfile).event_markers
    um = bioread.read_headers(ufile).event_markers
    assert um == cm


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


def clean_text(s):
    """
    Normalize line endings to \n and strip whitespace.

    These differences aren't important for our tests.
    """
    return s.replace('\r\n', '\n').strip()


class HTMLTextExtractor(HTMLParser):
    def handle_data(self, data):
        if not hasattr(self, 'content'):
            self.content = ''
        self.content += data



