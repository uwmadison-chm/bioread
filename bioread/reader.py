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
from contextlib import contextmanager

import numpy as np

import logging
# Re-adding the handler on reload causes duplicate log messages.
logger = logging.getLogger("bioread")
logger.setLevel(logging.WARNING)
log_handler = logging.StreamHandler()
log_handler.setLevel(logging.DEBUG)
log_handler.setFormatter(logging.Formatter("%(message)s"))
if len(logger.handlers) == 0:  # Avoid duplicate messages on reload
    logger.addHandler(log_handler)

import bioread.file_revisions as rev
from bioread import headers as bh
from bioread.headers import GraphHeader, ChannelHeader, ChannelDTypeHeader
from bioread.headers import ForeignHeader, MainCompressionHeader
from bioread.headers import ChannelCompressionHeader
from bioread.headers import PostMarkerHeader, V2JournalHeader, V4JournalHeader
from bioread.headers import V4JournalLengthHeader
from bioread.biopac import Datafile, EventMarker

CHUNK_SIZE = 1024 * 256  # A suggestion, probably not a terrible one.


class Reader(object):
    def __init__(self, acq_file=None):
        self.acq_file = acq_file
        self.datafile = None
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
        self.event_markers = None

    @classmethod
    def read(cls, fo, channel_indexes=None, target_chunk_size=CHUNK_SIZE):
        """ Read a biopac file into memory.

        fo: The name of the file to read, or a file-like object
        channel_indexes: The numbers of the channels you want to read
        target_chunk_size: The amount of data to read in a chunk.

        returns: reader.Reader.
        """
        with open_or_yield(fo, 'rb') as io:
            reader = cls(io)
            reader._read_headers()
            reader._read_data(channel_indexes, target_chunk_size)
        return reader

    @classmethod
    def read_headers(cls, fo):
        """ Read only the headers -- no data -- of a biopac file.
        """
        with open_or_yield(fo, 'rb') as io:
            reader = cls(io)
            reader._read_headers()
        return reader

    def stream(self, channel_indexes=None, target_chunk_size=CHUNK_SIZE):
        """ Set up and retun an iterator for streaming data.
        """
        if self.datafile is None:
            self._read_headers()
        if self.is_compressed:
            raise TypeError('Streaming is not supported for compressed files')
        self.acq_file.seek(self.data_start_offset)
        return make_chunk_reader(
            self.acq_file,
            self.datafile.channels,
            channel_indexes,
            target_chunk_size)

    @property
    def is_compressed(self):
        return self.graph_header.compressed

    def _read_headers(self):
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

        self.samples_per_second = 1000/self.graph_header.sample_time

        self.datafile = Datafile(
            graph_header=self.graph_header,
            channel_headers=self.channel_headers,
            foreign_header=self.foreign_header,
            channel_dtype_headers=self.channel_dtype_headers,
            samples_per_second=self.samples_per_second)

        self.data_length = self.datafile.data_length

        # In compressed files, markers come before compressed data. But
        # data_length is 0 for compressed files.
        self.marker_start_offset = (self.data_start_offset + self.data_length)
        self._read_markers()
        try:
            self._read_journal()
        except struct.error:
            logger.info("No journal information found.")
        if self.is_compressed:
            self.__read_compression_headers()

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
        self.journal = None
        self.journal_header = None
        if self.file_revision <= rev.V_400B:
            self.__read_journal_v2()
        else:
            self.__read_journal_v4()
        self.datafile.journal_header = self.journal_header
        self.datafile.journal = self.journal

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
            'utf-8').strip('\0')

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
            logger.debug("Reading {0} bytes of journal at {1}".format(
                self.journal_header.journal_len,
                self.acq_file.tell()))
            self.journal = self.acq_file.read(
                self.journal_header.journal_len).decode('utf-8').strip('\0')
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

    def _read_data(self, channel_indexes, target_chunk_size=CHUNK_SIZE):
        if self.is_compressed:
            self.__read_data_compressed(channel_indexes)
        else:
            self.__read_data_uncompressed(channel_indexes, target_chunk_size)

    def _read_markers(self):
        if self.marker_start_offset is None:
            self.read_headers()
        mh_class = bh.V2MarkerHeader
        mih_class = bh.V2MarkerItemHeader
        if self.file_revision >= rev.V_400B:
            mh_class = bh.V4MarkerHeader
            mih_class = bh.V4MarkerItemHeader
        self.marker_header = self.__single_header(
            self.marker_start_offset, mh_class)
        self.datafile.marker_header = self.marker_header
        self.__read_marker_items(mih_class)

    def __read_marker_items(self, marker_item_header_class):
        """
        self.acq_file must be seek()ed to the start of the first item header
        """
        event_markers = []
        marker_item_headers = []
        for i in range(self.marker_header.marker_count):
            mih = self.__single_header(
                self.acq_file.tell(), marker_item_header_class)
            marker_text_bytes = self.acq_file.read(mih.text_length)
            marker_text = marker_text_bytes.decode('utf-8').strip('\0')
            marker_item_headers.append(mih)
            marker_channel = self.datafile.channel_order_map.get(
                mih.channel_number)
            event_markers.append(EventMarker(
                sample_index=mih.sample_index,
                text=marker_text,
                channel_number=mih.channel_number,
                channel=marker_channel,
                type_code=mih.type_code))
        self.marker_item_headers = marker_item_headers
        self.datafile.marker_item_headers = marker_item_headers
        self.datafile.event_markers = event_markers

    def __read_data_compressed(self, channel_indexes):
        # At least in post-4.0 files, the compressed data isn't interleaved at
        # all. It's stored in uniform compressed blocks -- this probably
        # compresses far better than interleaved data.
        # Strangely, the compressed data seems to always be little-endian.
        if channel_indexes is None:
            channel_indexes = np.arange(len(self.datafile.channels))

        for i in channel_indexes:
            cch = self.channel_compression_headers[i]
            channel = self.datafile.channels[i]
            # Data seems to always be little-endian
            dt = channel.dtype.newbyteorder("<")
            self.acq_file.seek(cch.compressed_data_offset)
            comp_data = self.acq_file.read(cch.compressed_data_len)
            decomp_data = zlib.decompress(comp_data)
            channel.raw_data = np.fromstring(decomp_data, dtype=dt)

    def __read_data_uncompressed(self, channel_indexes, target_chunk_size):
        self.acq_file.seek(self.data_start_offset)
        # This will fill self.datafile.channels with data.
        read_uncompressed(
            self.acq_file,
            self.datafile.channels,
            channel_indexes,
            target_chunk_size)

    def __set_order_and_version(self):
        # Try unpacking the version string in both a bid and little-endian
        # fashion. Version string should be a small, positive integer.
        self.acq_file.seek(0)
        # No byte order flag -- we're gonna figure it out.
        gh = GraphHeader(rev.V_ALL, '')
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
        return "Reader('{0}')".format(self.acq_file)


@contextmanager
def open_or_yield(thing, mode):
    """ If 'thing' is a string, open it and yield it. Otherwise, yield it.

    This lets you use a filename, open file, other IO object. If 'thing' was
    a filename, the file is guaranteed to be closed after yielding.
    """
    if isinstance(thing, six.string_types):
        with open(thing, mode) as f:
            yield(f)
    else:
        yield(thing)


class ChunkBuffer(object):
    def __init__(self, channel):
        self.channel = channel
        self.buffer = None
        self.channel_slice = slice(0, 0)


def read_uncompressed(
        f,
        channels,
        channel_indexes=None,
        target_chunk_size=CHUNK_SIZE):
    """
    Read the uncompressed data.

    This function will read the data from an open IO object f (which must be
    seek()ed to the start of the data) into the raw_data attribute of
    channels.

    channel_indexes is a list of indexes of the channels we want to read
    (if None, read all the channels). Other channels' raw_data will be set
    to None.

    target_chunk_size gives a general idea of how much data the program should
    read into memory at a time. You can probably always leave this as at its
    default.

    This function returns nothing; it modifies channels in-place.

    Uncompressed data are stored in .acq files in an interleaved format --
    as the data streams off the amps, it's stored directly. So, with three
    channels, your data might look like (spaces added for clarity):

    012 012 012 012 012 ...

    Each channel can also have a frequency divider, which tells us this
    channel is recorded every nth occurence of the file's base sampling rate.

    If our three channels have frequency dividers [1, 4, 2], the data pattern
    would look like (again, with spaces between repetitions):
    0120020 0120020 0120020 ...

    """
    if channel_indexes is None:
        channel_indexes = np.arange(len(channels))

    for i in channel_indexes:
        channels[i]._allocate_raw_data()

    chunker = make_chunk_reader(
        f, channels, channel_indexes, target_chunk_size)
    for chunk_buffers in chunker:
        for i in channel_indexes:
            ch = channels[i]
            buf = chunk_buffers[i]
            logger.debug('Storing {0} samples to {1} of channel {2}'.format(
                len(buf.buffer), buf.channel_slice, i))
            ch.raw_data[buf.channel_slice] = buf.buffer[:]


def make_chunk_reader(
        f,
        channels,
        channel_indexes=None,
        target_chunk_size=CHUNK_SIZE):

    if channel_indexes is None:
        channel_indexes = np.arange(len(channels))

    byte_pattern = chunk_byte_pattern(channels, target_chunk_size)
    logger.debug('Using chunk size: {0} bytes'.format(len(byte_pattern)))
    buffers = [ChunkBuffer(c) for c in channels]
    return read_chunks(f, buffers, byte_pattern, channel_indexes)


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
        channel_bytes_remaining -= np.bincount(
            pat, minlength=len(channel_bytes_remaining))
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
