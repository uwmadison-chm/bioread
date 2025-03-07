# coding: utf8
# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2025 Board of Regents of the University of Wisconsin System
#
# Written Nate Vack <njvack@wisc.edu> with research from John Ollinger
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison
# Extended by Alexander Schlemmer.

from __future__ import with_statement, division
import struct
import zlib
from contextlib import contextmanager

import numpy as np

import bioread.file_revisions as rev
from bioread import headers as bh
from bioread.biopac import Datafile, EventMarker
from bioread.header_reader import HeaderReader
from bioread.dtype_header_reader import DTypeHeaderReader
from bioread.journal_reader import JournalReader
from bioread.marker_reader import MarkerReader
from bioread.data_reader import DataReader, CHUNK_SIZE

import logging
# Re-adding the handler on reload causes duplicate log messages.
logger = logging.getLogger("bioread")
logger.setLevel(logging.WARNING)
log_handler = logging.StreamHandler()
log_handler.setLevel(logging.DEBUG)
log_handler.setFormatter(logging.Formatter("%(message)s"))
if len(logger.handlers) == 0:  # Avoid duplicate messages on reload
    logger.addHandler(log_handler)


# How far past the foreign data header we're willing to go looking for the
# channel dtype headers
MAX_DTYPE_SCANS = 4096

READ_EXCEPTIONS = (
    ValueError,
    IOError,
    OSError,
    UnicodeDecodeError,
    EOFError,
    struct.error,
    IndexError
)


class Reader(object):
    def __init__(self, acq_file=None):
        self.acq_file = acq_file
        self.encoding = None  # We're gonna guess from _set_order_and_version
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
        self.data_length = None
        self.marker_start_offset = None
        self.marker_header = None
        self.marker_item_headers = None
        self.marker_metadata_pre_header = None
        self.marker_metadata_headers = None
        self.event_markers = None
        self.read_errors = []
        
        # Initialize readers
        self.header_reader = None
        self.dtype_header_reader = None
        self.journal_reader = None
        self.marker_reader = None
        self.data_reader = None

    # Public methods
    @classmethod
    def read(cls,
             fo,
             channel_indexes=None,
             target_chunk_size=CHUNK_SIZE):
        """ Read a biopac file into memory.

        fo: The name of the file to read, or a file-like object
        channel_indexes: The numbers of the channels you want to read
        target_chunk_size: The amount of data to read in a chunk.

        returns: reader.Reader.
        """
        with open_or_yield(fo, 'rb') as io:
            reader = cls(io)
            try:
                reader._read_headers()
            except READ_EXCEPTIONS as e:
                pass
            # We've already logged the error. Try to read the data.
            try:
                reader._read_data(channel_indexes, target_chunk_size)
            except READ_EXCEPTIONS as e:
                # Log and print the error, but consume the exception so the
                # caller still gets the reader.
                logger.error(f"Error reading data: {e}")
                reader.read_errors.append(str(e))
        return reader

    @classmethod
    def read_headers(cls, fo):
        """ Read only the headers -- no data -- of a biopac file.
        """
        with open_or_yield(fo, 'rb') as io:
            reader = cls(io)
            try:
                reader._read_headers()
            except READ_EXCEPTIONS as e:
                # Consume the exception so the caller still gets the reader.
                pass
        return reader

    def stream(self, channel_indexes=None, target_chunk_size=CHUNK_SIZE):
        """ Set up and retun an iterator for streaming data.
        """
        if self.datafile is None:
            self._read_headers()
        if self.is_compressed:
            raise TypeError('Streaming is not supported for compressed files')
            
        # Initialize data reader if needed
        if self.data_reader is None:
            self.data_reader = DataReader(
                self.acq_file, 
                self.datafile, 
                self.data_start_offset
            )
            
        return self.data_reader.stream(channel_indexes, target_chunk_size)

    @property
    def is_compressed(self):
        return self.graph_header.compressed

    def __repr__(self):
        return "Reader('{0}')".format(self.acq_file)

    # Header reading methods
    def _read_headers(self):
        logger.debug("I am in _read_headers")
        
        # Initialize header reader
        self.header_reader = HeaderReader(self.acq_file, self.byte_order_char, 
                                         self.file_revision, self.encoding)
        
        # Set byte order and version if not already set
        if self.byte_order_char is None:
            self.file_revision, self.byte_order_char, self.encoding = self.header_reader.set_order_and_version()

        # Initialize specialized readers
        self.dtype_header_reader = DTypeHeaderReader(self.header_reader)
        self.journal_reader = JournalReader.create_journal_reader(self.header_reader, self.file_revision)
        self.marker_reader = MarkerReader.create_marker_reader(self.header_reader, self.file_revision)

        # Read graph header
        graph_header_class = bh.get_graph_header_class(self.file_revision)
        self.graph_header = self.header_reader.single_header(0, graph_header_class)
        channel_count = self.graph_header.channel_count

        # Read padding headers
        pad_start = self.graph_header.effective_len_bytes
        pad_headers = self.header_reader.multi_headers(
            self.graph_header.expected_padding_headers,
            pad_start,
            bh.UnknownPaddingHeader)
        
        # Read channel headers
        ch_start = pad_start + sum(
            [ph.effective_len_bytes for ph in pad_headers])
        # skip past the unknown padding header
        _ = self.header_reader.single_header(ch_start, bh.UnknownPaddingHeader)
        channel_header_class = bh.get_channel_header_class(self.file_revision)
        self.channel_headers = self.header_reader.multi_headers(channel_count,
                                                  ch_start, channel_header_class)
        ch_len = self.channel_headers[0].effective_len_bytes

        for i, ch in enumerate(self.channel_headers):
            logger.debug("Channel header %s: %s" % (i, ch.data))

        # Read foreign header
        fh_start = ch_start + len(self.channel_headers)*ch_len
        foreign_header_class = bh.get_foreign_header_class(self.file_revision)
        self.foreign_header = self.header_reader.single_header(fh_start, foreign_header_class)

        # Read channel dtype headers
        cdh_start = fh_start + self.foreign_header.effective_len_bytes
        self.channel_dtype_headers, self.data_start_offset = self.dtype_header_reader.scan_for_dtype_headers(
            cdh_start, channel_count)
        if self.channel_dtype_headers is None:
            raise ValueError("Can't find valid channel data type headers")

        for i, cdt in enumerate(self.channel_dtype_headers):
            logger.debug("Channel %s: type_code: %s, offset: %s" % (
                i, cdt.type_code, cdt.offset
            ))

        logger.debug("Computed data start offset: %s" % self.data_start_offset)

        # Calculate samples per second
        self.samples_per_second = 1000/self.graph_header.sample_time

        # Create datafile
        logger.debug("About to allocate a Datafile")
        self.datafile = Datafile(
            graph_header=self.graph_header,
            channel_headers=self.channel_headers,
            foreign_header=self.foreign_header,
            channel_dtype_headers=self.channel_dtype_headers,
            samples_per_second=self.samples_per_second)

        logger.debug("Allocated a datafile!")

        # Calculate data length
        self.data_length = self.datafile.data_length
        logger.debug("Computed data length: %s" % self.data_length)

        # In compressed files, markers come before compressed data. But
        # data_length is 0 for compressed files.
        self.marker_start_offset = (self.data_start_offset + self.data_length)
        
        # Read markers
        self._read_markers()
        
        # Read journal
        try:
            self._read_journal()
        except struct.error:
            logger.info("No journal information found.")
            
        # Read compression headers if needed
        if self.is_compressed:
            self._read_compression_headers()
            
        # Initialize data reader
        self.data_reader = DataReader(
            self.acq_file, 
            self.datafile, 
            self.data_start_offset
        )
        
        # Store compression headers in datafile
        if self.is_compressed:
            self.datafile.main_compression_header = self.main_compression_header
            self.datafile.channel_compression_headers = self.channel_compression_headers

    # Journal reading methods
    def _read_journal(self):
        self.journal = None
        self.journal_header = None

        try:
            if self.file_revision <= rev.V_400B:
                self.journal_header, self.journal = self.journal_reader.read_journal()
            else:
                self.journal_header, self.journal, self.journal_length_header = self.journal_reader.read_journal()
        except READ_EXCEPTIONS as e:
            logger.error(f"Error reading journal: {e}")
            self.read_errors.append(f"Error reading journal: {e}")
            raise
            
        self.datafile.journal_header = self.journal_header
        self.datafile.journal = self.journal

    # Marker reading methods
    def _read_markers(self):
        if self.marker_start_offset is None:
            self.read_headers()
            
        # Read markers using the appropriate marker reader
        result = self.marker_reader.read_markers(self.marker_start_offset, self.graph_header)
        
        # Unpack the result
        self.marker_header, self.marker_item_headers, event_markers, self.marker_metadata_pre_header, self.marker_metadata_headers = result
            
        # Update channel references in event markers
        for marker in event_markers:
            marker.channel = self.datafile.channel_order_map.get(
                marker.channel_number)
                
        self.datafile.marker_header = self.marker_header
        self.datafile.marker_item_headers = self.marker_item_headers
        self.datafile.event_markers = event_markers
        self.event_markers = event_markers

    def _read_compression_headers(self):
        # We need to have read the markers and journal; this puts us
        # at the correct file offset.
        self.marker_start_offset = self.data_start_offset
        main_ch_start = self.acq_file.tell()
        main_compression_header_class = bh.get_main_compression_header_class(self.file_revision)
        self.main_compression_header = self.header_reader.single_header(
            main_ch_start, main_compression_header_class)
        cch_start = (main_ch_start +
                     self.main_compression_header.effective_len_bytes)
        self.channel_compression_headers = self.header_reader.multi_headers(
            self.graph_header.channel_count, cch_start,
            bh.ChannelCompressionHeader)

    def _read_data(self, channel_indexes, target_chunk_size=CHUNK_SIZE):
        if self.data_reader is None:
            self.data_reader = DataReader(
                self.acq_file, 
                self.datafile, 
                self.data_start_offset
            )
            
        self.data_reader.read_data(channel_indexes, target_chunk_size)


@contextmanager
def open_or_yield(thing, mode):
    """ If 'thing' is a string, open it and yield it. Otherwise, yield it.

    This lets you use a filename, open file, other IO object. If 'thing' was
    a filename, the file is guaranteed to be closed after yielding.
    """
    if isinstance(thing, str):
        logger.debug(f"Opening file: {thing}")
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
        chunk_data = np.frombuffer(
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
    pattern_mask = np.zeros(len(byte_pattern), dtype=bool)
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
