# coding: utf8
# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2025 Board of Regents of the University of Wisconsin System
#
# Written Nate Vack <njvack@wisc.edu> with research from John Ollinger
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison
# Extended by Alexander Schlemmer.

import logging
import zlib
import numpy as np

logger = logging.getLogger("bioread")

# This is how much interleaved uncompressed data we'll read at a time.
CHUNK_SIZE = 1024 * 256  # A suggestion, probably not a terrible one.


class ChunkBuffer:
    def __init__(self, channel):
        self.channel = channel
        self.buffer = None
        self.channel_slice = slice(0, 0)


class DataReader:
    """
    A class to handle reading data from a BIOPAC file.
    """

    def __init__(self, acq_file, datafile, data_start_offset):
        """
        Initialize a DataReader.

        Parameters
        ----------
        acq_file : file-like object
            The file to read data from
        datafile : Datafile
            The datafile object to store data in
        data_start_offset : int
            The offset where the data starts
        """
        self.acq_file = acq_file
        self.datafile = datafile
        self.data_start_offset = data_start_offset

    def read_data(self, channel_indexes=None, target_chunk_size=CHUNK_SIZE):
        """
        Read data from the file.

        Parameters
        ----------
        channel_indexes : list, optional
            The indexes of the channels to read
        target_chunk_size : int, optional
            The target chunk size for reading uncompressed data
        """
        if self.datafile.is_compressed:
            self.read_data_compressed(channel_indexes)
        else:
            self.read_data_uncompressed(channel_indexes, target_chunk_size)

    def read_data_compressed(self, channel_indexes=None):
        """
        Read compressed data from the file.

        Parameters
        ----------
        channel_indexes : list, optional
            The indexes of the channels to read
        """
        # The compressed data isn't interleaved at
        # all. It's stored in uniform compressed blocks -- this probably
        # compresses far better than interleaved data.
        # Strangely, the compressed data seems to always be little-endian.
        if channel_indexes is None:
            channel_indexes = np.arange(len(self.datafile.channels))

        for i in channel_indexes:
            cch = self.datafile.channel_compression_headers[i]
            channel = self.datafile.channels[i]
            self.acq_file.seek(cch.compressed_data_offset)
            comp_data = self.acq_file.read(cch.compressed_data_len)
            decomp_data = zlib.decompress(comp_data)
            channel.raw_data = np.frombuffer(decomp_data, dtype=channel.dtype)

    def read_data_uncompressed(
        self, channel_indexes=None, target_chunk_size=CHUNK_SIZE
    ):
        """
        Read uncompressed data from the file.

        Parameters
        ----------
        channel_indexes : list, optional
            The indexes of the channels to read
        target_chunk_size : int, optional
            The target chunk size for reading
        """
        self.acq_file.seek(self.data_start_offset)
        # This will fill self.datafile.channels with data.
        read_uncompressed(
            self.acq_file, self.datafile.channels, channel_indexes, target_chunk_size
        )

    def stream(self, channel_indexes=None, target_chunk_size=CHUNK_SIZE):
        """
        Set up a streaming iterator for the data.

        Parameters
        ----------
        channel_indexes : list, optional
            The indexes of the channels to read
        target_chunk_size : int, optional
            The target chunk size for reading

        Returns
        -------
        iterator
            An iterator that yields chunks of data
        """
        if self.datafile.is_compressed:
            raise TypeError("Streaming is not supported for compressed files")
        self.acq_file.seek(self.data_start_offset)
        return make_chunk_reader(
            self.acq_file, self.datafile.channels, channel_indexes, target_chunk_size
        )


# Utility functions for data reading


def read_uncompressed(f, channels, channel_indexes=None, target_chunk_size=CHUNK_SIZE):
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

    chunker = make_chunk_reader(f, channels, channel_indexes, target_chunk_size)
    for chunk_buffers in chunker:
        for i in channel_indexes:
            ch = channels[i]
            buf = chunk_buffers[i]
            logger.debug(
                "Storing {0} samples to {1} of channel {2}".format(
                    len(buf.buffer), buf.channel_slice, i
                )
            )
            ch.raw_data[buf.channel_slice] = buf.buffer[:]


def make_chunk_reader(f, channels, channel_indexes=None, target_chunk_size=CHUNK_SIZE):
    if channel_indexes is None:
        channel_indexes = np.arange(len(channels))

    byte_pattern = chunk_byte_pattern(channels, target_chunk_size)
    logger.debug("Using chunk size: {0} bytes".format(len(byte_pattern)))
    buffers = [ChunkBuffer(c) for c in channels]
    return read_chunks(f, buffers, byte_pattern, channel_indexes)


def read_chunks(f, buffers, byte_pattern, channel_indexes):
    """
    Read data in chunks from f. For each chunk, yield a list of buffers with
    information on how much of the buffer is filled and where the data should
    go in the target array.
    """
    channel_bytes_remaining = np.array([b.channel.data_length for b in buffers])
    chunk_number = 0
    while np.sum(channel_bytes_remaining) > 0:
        pat = chunk_pattern(byte_pattern, channel_bytes_remaining)
        chunk_bytes = len(pat)
        logger.debug(
            "Chunk {0}: {1} bytes at {2}".format(chunk_number, chunk_bytes, f.tell())
        )
        chunk_data = np.frombuffer(f.read(chunk_bytes), dtype="b", count=chunk_bytes)
        update_buffers_with_data(chunk_data, buffers, pat, channel_indexes)

        yield buffers
        channel_bytes_remaining -= np.bincount(
            pat, minlength=len(channel_bytes_remaining)
        )
        logger.debug("Channel bytes remaining: {0}".format(channel_bytes_remaining))
        chunk_number += 1


def chunk_pattern(byte_pattern, channel_bytes_remaining):
    """Trim a byte pattern depending on how many bytes remain in each channel.

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
    trimmed_pattern = byte_pattern[0 : len(data)]
    for i in channel_indexes:
        buf = buffers[i]
        buf.buffer = data[trimmed_pattern == i]
        buf.buffer.dtype = buf.channel.dtype
        old_slice = buf.channel_slice
        buf.channel_slice = slice(old_slice.stop, old_slice.stop + len(buf.buffer))


def chunk_byte_pattern(channels, target_chunk_size):
    """Compute a byte layout for a chunk of data.

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
    """Compute the pattern of samples in a file's uncompressed data stream.

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
    pattern_slots = (
        np.arange(base_len).repeat(channel_count).reshape((base_len, channel_count))
    )
    pattern_mask = (pattern_slots % dividers) == 0
    channel_slots = np.tile(np.arange(channel_count), (base_len, 1))
    return channel_slots[pattern_mask]


def chunk_pattern_reps(target_chunk_size, pattern_byte_length):
    """
    The number of times we'll actually repeat the pattern in a chunk.
    Must always be at least 1.
    """
    return max(1, target_chunk_size // pattern_byte_length)


def least_common_multiple(*ar):
    """Compute least common multiple of n numbers.

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
    """Iterative method to compute greatest common denominator."""
    while not b == 0:
        a, b = b, a % b
    return a
