# coding: utf8
# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2016 Board of Regents of the University of Wisconsin System
#
# Written Nate Vack <njvack@wisc.edu> with research from John Ollinger
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison
# Project home: http://github.com/njvack/bioread
# Extended by Alexander Schlemmer.

from __future__ import with_statement, division
from bioread.vendor import six
import struct
import zlib

import numpy as np

import logging
# Re-adding the handler on reload causes duplicate log messages.
logger = logging.getLogger("bioread")
logger.setLevel(logging.DEBUG)
log_handler = logging.StreamHandler()
log_handler.setLevel(logging.DEBUG)
log_handler.setFormatter(logging.Formatter("%(message)s"))
if len(logger.handlers) == 0:  # Avoid duplicate messages on reload
    logger.addHandler(log_handler)

from bioread.file_revisions import *
from bioread import headers as bh
from bioread.headers import GraphHeader, ChannelHeader, ChannelDTypeHeader
from bioread.headers import ForeignHeader, MainCompressionHeader
from bioread.headers import ChannelCompressionHeader
from bioread.headers import PostMarkerHeader, V2JournalHeader, V4JournalHeader
from bioread.headers import V4JournalLengthHeader
from bioread.biopac import Datafile, Channel, Marker


class AcqReader(object):
    """
    Main class for reading AcqKnowledge files. You'll probably call it like:
    >>> data = AcqReader.read("some_file.acq")
    """

    def __init__(self, acq_file):
        self.acq_file = acq_file
        # This must be set by _set_order_and_version
        self.byte_order_char = None
        self.file_revision = None
        self.samples_per_second = None
        self.graph_header = None
        self.channel_headers = []
        self.foreign_header = None
        self.channel_dtype_headers = []
        self.main_compression_header = None
        self.channel_compression_headers = []
        self.data_start_offset = None
        self.marker_start_offset = None
        self.marker_header = None
        self.marker_item_headers = None
        self.markers = None

    @classmethod
    def read_file(cls, fo):
        """
        The main method to quickly read a biopac file into memory.

        fo: The name of the file to read, or a file-like object

        returns: biopac.Datafile
        """
        if isinstance(fo, six.string_types):
            with open(fo, 'rb') as f:
                reader = cls(f)
                return reader.read()
        else:
            reader = cls(fo)
            return reader.read()

    @classmethod
    def read_without_data(cls, fo):
        if isinstance(fo, six.string_types):
            with open(fo, 'rb') as f:
                reader = cls(f)
                return reader.read_headers()
        else:
            reader = cls(fo)
            return reader.read_headers()

    def read(self):
        df = self.read_headers()

        self._read_data()
        df.channels = self.channels
        df.marker_header = self.marker_header
        df.marker_item_headers = self.marker_item_headers
        df.markers = self.markers
        self.data_file = df
        return self.data_file

    @property
    def is_compressed(self):
        return self.graph_header.compressed

    def read_headers(self):
        if self.byte_order_char is None:
            self.__set_order_and_version()

        self.graph_header = self.__single_header(0, GraphHeader)
        channel_count = self.graph_header.channel_count

        ch_start = self.graph_header.effective_len_bytes
        self.channel_headers = self.__multi_headers(channel_count,
                                                    ch_start, ChannelHeader)
        ch_len = self.channel_headers[0].effective_len_bytes

        fh_start = ch_start + len(self.channel_headers)*ch_len
        self.foreign_header = self.__single_header(fh_start, ForeignHeader)

        cdh_start = fh_start + self.foreign_header.effective_len_bytes
        self.channel_dtype_headers = self.__multi_headers(
            channel_count, cdh_start, ChannelDTypeHeader)
        cdh_len = self.channel_dtype_headers[0].effective_len_bytes

        self.data_start_offset = (cdh_start + (cdh_len * channel_count))
        data_length = sum(
            [c.point_count * cd.sample_size for c, cd in
                zip(self.channel_headers, self.channel_dtype_headers)])
        self.marker_start_offset = self.data_start_offset + data_length
        if self.is_compressed:
            # In this case, the marker and journal come *before* the data
            self.marker_start_offset = self.data_start_offset
        self._read_markers()
        self._read_journal()
        if self.is_compressed:
            self.__read_compression_headers()

        self.samples_per_second = 1000/self.graph_header.sample_time
        return Datafile(
            graph_header=self.graph_header,
            channel_headers=self.channel_headers,
            foreign_header=self.foreign_header,
            channel_dtype_headers=self.channel_dtype_headers,
            samples_per_second=self.samples_per_second)

    def __read_compression_headers(self):
        # We need to have read the markers and journal; this puts us
        # at the correct file offset.
        self.marker_start_offset = self.data_start_offset
        main_ch_start = self.acq_file.tell()
        self.main_compression_header = self.__single_header(
            main_ch_start, MainCompressionHeader)
        cch_start = (main_ch_start +
                     self.main_compression_header.effective_len_bytes)
        self.channel_compression_headers = self.__multi_headers(
            self.graph_header.channel_count, cch_start,
            ChannelCompressionHeader)

    def _read_journal(self):
        if self.file_revision <= V_400B:
            self.__read_journal_v2()
        else:
            self.__read_journal_v4()

    def __read_journal_v2(self):
        self.post_marker_header = self.__single_header(
            self.acq_file.tell(), PostMarkerHeader)
        logger.debug(self.acq_file.tell())
        logger.debug(self.post_marker_header.rep_bytes)
        self.acq_file.seek(self.post_marker_header.rep_bytes, 1)
        logger.debug(self.acq_file.tell())
        self.journal_header = self.__single_header(
            self.acq_file.tell(), V2JournalHeader)
        self.journal = self.acq_file.read(
            self.journal_header.data['lJournalLen']).decode(
            'utf-8').strip('\x00')

    def __read_journal_v4(self):
        self.journal_length_header = self.__single_header(
            self.acq_file.tell(),
            V4JournalLengthHeader)
        journal_len = self.journal_length_header.journal_len
        self.journal = None
        jh = V4JournalHeader(
            self.file_revision, self.byte_order_char)
        # If journal_length_header.journal_len is small, we don't have a
        # journal to read.
        if (jh.effective_len_bytes <= journal_len):
            self.journal_header = self.__single_header(
                self.acq_file.tell(),
                V4JournalHeader)
            logger.debug("Reading {0} bytes of journal at {0}".format(
                self.journal_header.journal_len,
                self.acq_file.tell()))
            self.journal = self.acq_file.read(
                self.journal_header.journal_len).decode('utf-8').strip('\x00')
        # Either way, we should seek to this point.
        self.acq_file.seek(self.journal_length_header.data_end)

    def __single_header(self, start_offset, h_class):
        return self.__multi_headers(1, start_offset, h_class)[0]

    def __multi_headers(self, num, start_offset, h_class):
        headers = []
        last_h_len = 0  # This will be changed reading the channel headers
        h_offset = start_offset
        for i in range(num):
            h_offset += last_h_len
            logger.debug(
                "Reading {0} at offset {1}".format(h_class, h_offset))
            h = h_class(self.file_revision, self.byte_order_char)
            h.unpack_from_file(self.acq_file, h_offset)
            last_h_len = h.effective_len_bytes
            headers.append(h)
        return headers

    def __build_channels(self):
        # Build empty channels, ready to get data from the file.
        return [
            Channel.from_headers(ch, cdh, self.samples_per_second)
            for ch, cdh in
            zip(self.channel_headers, self.channel_dtype_headers)
        ]

    def _read_data(self):
        self.channels = self.__build_channels()
        if self.is_compressed:
            self.__read_data_compressed(self.channels)
        else:
            self.__read_data_uncompressed(self.channels)

    def _read_markers(self):
        if self.marker_start_offset is None:
            self.read_headers()
        mh_class = bh.V2MarkerHeader
        mih_class = bh.V2MarkerItemHeader
        if self.file_revision >= V_400B:
            mh_class = bh.V4MarkerHeader
            mih_class = bh.V4MarkerItemHeader
        self.marker_header = self.__single_header(
            self.marker_start_offset, mh_class)
        self.__read_marker_items(mih_class)

    def __read_marker_items(self, marker_item_header_class):
        """
        self.acq_file must be seek()ed to the start of the first item header
        """
        self.markers = []
        self.marker_item_headers = []
        for i in range(self.marker_header.marker_count):
            mih = self.__single_header(
                self.acq_file.tell(), marker_item_header_class)
            marker_text_bytes = self.acq_file.read(mih.text_length)
            marker_text = marker_text_bytes.decode('utf-8').strip('\0')
            self.marker_item_headers.append(mih)
            self.markers.append(Marker(
                mih.sample_index, marker_text, mih.channel, mih.style))

    def __read_data_compressed(self, channels):
        # At least in post-4.0 files, the compressed data isn't interleaved at
        # all. It's stored in uniform compressed blocks -- this probably
        # compresses far better than interleaved data.
        # Strangely, the compressed data seems to always be little-endian.
        for i in range(len(channels)):
            cch = self.channel_compression_headers[i]
            chan = channels[i]
            # Data seems to always be little-endian
            dt = chan.dtype.newbyteorder("<")
            self.acq_file.seek(cch.compressed_data_offset)
            comp_data = self.acq_file.read(cch.compressed_data_len)
            decomp_data = zlib.decompress(comp_data)
            chan.raw_data = np.fromstring(decomp_data, dtype=dt)

    def __read_data_uncompressed(self, channels):
        # The data in the file are interleaved, so we'll potentially have
        # a different amount of data to read at each time slice.
        # It's possible we won't have any data for some time slices, I think.
        # The BIOPAC engineers tell you not to even try reading interleaved
        # data. Wusses.

        # Using adapted algorithm by Sven Marnarch from:
        # http://stackoverflow.com/questions/4227990
        self.stream_sample_indexes = self.__stream_sample_indexes(channels)
        self.samples_per_block = len(self.stream_sample_indexes)
        self.total_samples = sum([c.point_count for c in channels])
        self.total_blocks = int(
            np.ceil(float(self.total_samples)/self.samples_per_block))

        self.all_sample_indexes = np.tile(self.stream_sample_indexes,
            self.total_blocks)
        self.channel_lengths = np.array([c.point_count for c in channels])
        self.channel_sizes = np.array([c.sample_size for c in channels])
        self.sample_counts = self.__sample_counts(
            self.stream_sample_indexes, self.total_blocks)
        self.sample_mask = (
            self.sample_counts < self.channel_lengths[self.all_sample_indexes])
        self.sample_map = self.all_sample_indexes[self.sample_mask]
        # The mapping of actual bytes on disk to channels.
        self.data_map = self.sample_map.repeat(
            self.channel_sizes[self.sample_map])

        self.acq_file.seek(self.data_start_offset)
        self.buf = np.fromfile(self.acq_file, np.ubyte, len(self.data_map))
        for i, ch in enumerate(channels):
            self.__copy_uncompressed_data(ch, i, self.data_map, self.buf)

    def __stream_sample_indexes(self, channels):
        """
        Returns the shortest repeating pattern of samples that'll appear in
        our data stream. If our frequency_dividers look like [1,2,4], we'll return
        [0,1,2,0,0,1,0]
        """
        dividers = [c.frequency_divider for c in channels]
        channel_lcm = least_common_multiple(*dividers)
        # Make a list like [0,1,2,0,0,1,0]
        stream_sample_indexes = [
            ch_idx for pat_idx in range(channel_lcm)
            for ch_idx, div in enumerate(dividers)
            if pat_idx % div == 0]
        return np.array(stream_sample_indexes, dtype=np.int32)

    def __copy_uncompressed_data(self, channel, d_map, channel_index, buf):
        mask = d_map == channel_index
        logger.debug(mask)
        logger.debug(len(mask))
        logger.debug(len(buf))
        channel.raw_data = buf[mask]
        channel.raw_data.dtype = channel.dtype

    def __to_byte_indexes(self, sample_indexes, channels):
        """
        Transform an array of sample_indexes (eg [0,1,2,0,0,1,0]) into an
        array of byte indexes: eg [0,0,1,1,2,2,2,2,2,2,2,2,0,0,0,0,1,1,0,0]
        if our channels' sample sizes are [2,2,8]
        """
        repeats = [channels[i].sample_size for i in sample_indexes]
        return np.array(sample_indexes).repeat(repeats)

    def __sample_counts(self, sample_indexes, reps):
        tile_bins = np.histogram(sample_indexes, np.max(sample_indexes)+1)[0]
        tile_mult = tile_bins[sample_indexes]
        first_steps = self.__running_counts(sample_indexes)
        tiled = np.tile(tile_mult, reps).reshape(reps, -1)
        multiplier = np.reshape(np.arange(reps, dtype=np.int32), (reps, 1))
        tiled *= multiplier
        tiled += first_steps
        return tiled.ravel()

    def __running_counts(self, sample_indexes):
        uniques = np.unique(sample_indexes)
        my_range = np.arange(len(sample_indexes), dtype=np.int32)
        counts = np.empty(sample_indexes.shape, dtype=np.int32)
        for i in uniques:
            counts[sample_indexes == i] = my_range[i]
        return counts

    def __set_order_and_version(self):
        # Try unpacking the version string in both a bid and little-endian
        # fashion. Version string should be a small, positive integer.
        self.acq_file.seek(0)
        # No byte order flag -- we're gonna figure it out.
        gh = GraphHeader(V_ALL, '')
        ver_fmt_str = gh.format_string
        ver_len = struct.calcsize('<'+ver_fmt_str)
        ver_data = self.acq_file.read(ver_len)

        byte_order_chars = ['<', '>']
        # Try both ways.
        byte_order_versions = [
            (struct.unpack(boc+ver_fmt_str, ver_data)[1], boc)
            for boc in byte_order_chars
        ]

        # Limit to positive numbers, choose smallest.
        byte_order_versions = sorted([
            bp for bp in byte_order_versions if bp[0] > 0])
        bp = byte_order_versions[0]

        self.byte_order_char = bp[1]
        self.file_revision = bp[0]

    def __repr__(self):
        return "AcqReader('{0}')".format(self.acq_file)


class ChunkBuffer(object):
    def __init__(self, channel):
        self.channel = channel
        self.buffer = None
        self.channel_slice = slice(0, 0)


def read_uncompressed(
        f,
        target_chunk_size,
        channels,
        channel_indexes=None):
    """
    Read the uncompressed data.
    """
    channel_indexes = channel_indexes or np.arange(len(channels))
    for i in channel_indexes:
        channels[i]._allocate_raw_data()
    byte_pattern = chunk_byte_pattern(channels, target_chunk_size)
    logger.debug('Using chunk size: {0} bytes'.format(len(byte_pattern)))
    buffers = [ChunkBuffer(c) for c in channels]

    chunker = read_chunks(f, buffers, byte_pattern, channel_indexes)
    for chunk_buffers in chunker:
        for i in channel_indexes:
            ch = channels[i]
            buf = chunk_buffers[i]
            logger.debug('Storing {0} samples to {1} of channel {2}'.format(
                len(buf.buffer), buf.channel_slice, i))
            ch.raw_data[buf.channel_slice] = buf.buffer[:]


def read_chunks(f, buffers, byte_pattern, channel_indexes):
    """
    Read data in chunks from f. For each chunk, yield a list of buffers with
    information on how much of the buffer is filled and where the data should
    go in the target array.
    """
    channel_bytes_remaining = np.array(
        [b.channel.data_length for b in buffers])
    chunk_number = 0
    while np.sum(channel_bytes_remaining) > 0:
        pat = chunk_pattern(byte_pattern, channel_bytes_remaining)
        chunk_bytes = len(pat)
        logger.debug('Chunk {0}: {1} bytes at {2}'.format(
            chunk_number, chunk_bytes, f.tell()))
        chunk_data = np.fromstring(
            f.read(chunk_bytes), dtype="b", count=chunk_bytes)
        update_buffers_with_data(
            chunk_data, buffers, pat, channel_indexes)

        yield buffers
        channel_bytes_remaining -= np.bincount(pat)
        logger.debug('Channel bytes remaining: {0}'.format(
            channel_bytes_remaining))
        chunk_number += 1


def chunk_pattern(byte_pattern, channel_bytes_remaining):
    """ Trim a byte pattern depending on how many bytes remain in each channel.

    For some reason, the data at the end of the file doesn't work like you'd
    expect. You can, for example, be missing an expected sample in a slow-
    sampling channel.

    The solution is to use the number of bytes in a channel to determine the
    actual layout of the chunk.
    """
    # This is the normal case, we don't need to do anything.
    if np.all(np.bincount(byte_pattern) <= channel_bytes_remaining):
        return byte_pattern
    # For each channel, compute a set of indexes where we expect data.
    channel_byte_indexes = [
        np.where(byte_pattern == i)[0][0:rem]
        for i, rem in enumerate(channel_bytes_remaining)
    ]
    all_byte_indexes = np.concatenate(channel_byte_indexes)
    pattern_mask = np.zeros(len(byte_pattern), dtype=np.bool)
    pattern_mask[all_byte_indexes] = True
    return byte_pattern[pattern_mask]


def update_buffers_with_data(data, buffers, byte_pattern, channel_indexes):
    """
    Updates buffers with information from data. Returns nothing, modifies
    buffers in-place.
    """
    trimmed_pattern = byte_pattern[0:len(data)]
    for i in channel_indexes:
        buf = buffers[i]
        buf.buffer = data[trimmed_pattern == i]
        buf.buffer.dtype = buf.channel.dtype
        old_slice = buf.channel_slice
        buf.channel_slice = slice(
            old_slice.stop, old_slice.stop + len(buf.buffer))


def chunk_byte_pattern(channels, target_chunk_size):
    """ Compute a byte layout for a chunk of data.

    This pattern is the main thing we actually need -- from it, we can know
    how to make individual buffers and how much data to read.

    The actual chunk size will always be a multiple of the byte pattern
    length, and will generally be very close to target_chunk_size. Usually, it
    will be larger.
    """
    divs = np.array([c.frequency_divider for c in channels])
    sizes = np.array([c.sample_size for c in channels])
    spat = sample_pattern(divs)
    byte_counts = sizes[spat]  # Returns array the length of spat
    bpat = spat.repeat(byte_counts)
    reps = chunk_pattern_reps(target_chunk_size, len(bpat))
    return np.tile(bpat, reps)


def sample_pattern(frequency_dividers):
    """ Compute the pattern of samples in a file's uncompressed data stream.

    The basic algorithm:
    * Take the least common multiple of the frequency dividers. This is the
      "base" of the pattern length -- the most times a channel could appear in
      the pattern.
    * Make a [base_len x num_channels] dimension matrix, counting from 0 to
      pattern_len in each row -- call this "pattern_slots"
    * Make a pattern_mask -- a boolean mask where each channel slots modulo
      frequency_divider == 0
    * The pattern, then, are the pattern_slots where pattern_mask is true

    Note that this is not quite the byte pattern -- these samples can either
    be int16 or float64.
    """
    dividers = np.array(frequency_dividers)
    channel_count = len(dividers)
    base_len = least_common_multiple(*dividers)
    pattern_slots = np.arange(
        base_len).repeat(
        channel_count).reshape(
        (base_len, channel_count))
    pattern_mask = ((pattern_slots % dividers) == 0)
    channel_slots = np.tile(np.arange(channel_count), (base_len, 1))
    return channel_slots[pattern_mask]


def chunk_pattern_reps(target_chunk_size, pattern_byte_length):
    """
    The number of times we'll actually repeat the pattern in a chunk.
    Must always be at least 1.
    """
    return max(1, target_chunk_size // pattern_byte_length)


def least_common_multiple(*ar):
    """ Compute least common multiple of n numbers.

    Adapted from:
    http://stackoverflow.com/questions/147515/least-common-multiple-for-3-or-more-numbers

    Used in computing the repeating pattern of multichannel data that's
    sampled at different rates in each channel.
    """

    if len(ar) > 2:
        return least_common_multiple(ar[0], least_common_multiple(*ar[1:]))
    elif len(ar) == 2:
        return (ar[0] * ar[1]) // greatest_common_denominator(ar[0], ar[1])
    else:
        return ar[0]


def greatest_common_denominator(a, b):
    """ Iterative method to compute greatest common denominator. """
    while not b == 0:
        a, b = b, a % b
    return a
