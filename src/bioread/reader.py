# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2025 Board of Regents of the University of Wisconsin System
#
# Written Nate Vack <njvack@wisc.edu> with research from John Ollinger
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison
# Extended by Alexander Schlemmer.

import logging
import os
import struct
from contextlib import contextmanager
from io import IOBase
from typing import BinaryIO

from bioread import headers as bh
from bioread.biopac import Datafile
from bioread.data_reader import CHUNK_SIZE, DataReader
from bioread.file_revisions import version_string_guess
from bioread.header_reader import (
    ChannelDTypeHeaderReader,
    GraphHeaderReader,
    HeaderReader,
)
from bioread.journal_reader import JournalReader
from bioread.marker_reader import MarkerReader

logger = logging.getLogger(__name__)


READ_EXCEPTIONS = (
    ValueError,
    IOError,
    OSError,
    UnicodeDecodeError,
    EOFError,
    struct.error,
    IndexError,
)


def _describe_error(what, e):
    """Describe a read error, including any notes on where it happened."""
    notes = getattr(e, "__notes__", [])
    return "; ".join([f"Error reading {what}: {e}", *notes])


class Reader:
    def __init__(self, acq_file: BinaryIO):
        self.acq_file = acq_file
        self.encoding = None  # We're gonna guess from _set_order_and_version
        self.datafile: Datafile | None = None
        # This must be set by _set_order_and_version
        self.byte_order_char = None
        self.file_revision = None
        self.version_string = None
        self.samples_per_second = None
        self.headers = []
        self.graph_header: bh.GraphHeader | None = None
        self.channel_headers = []
        self.channel_dtype_headers = []
        self.channel_compression_headers = []
        self.data_start_offset = None
        self.data_length = None
        self.event_markers = None
        self.read_errors = []

    # Public methods
    @classmethod
    def read(
        cls,
        fo: str | os.PathLike | BinaryIO,
        channel_indexes: list[int] | None = None,
        target_chunk_size: int = CHUNK_SIZE,
    ) -> "Reader":
        """Read a biopac file into memory.

        fo: The name of the file to read, or a file-like object
        channel_indexes: The numbers of the channels you want to read
        target_chunk_size: The amount of data to read in a chunk.

        returns: reader.Reader.
        """
        with open_or_yield(fo, "rb") as io:
            reader = cls(io)
            try:
                reader._read_headers()
            except READ_EXCEPTIONS as e:
                # Keep going: often enough was read to get at the data.
                reader._record_read_error("headers", e)
            try:
                reader._read_data(channel_indexes, target_chunk_size)
            except READ_EXCEPTIONS as e:
                reader.read_errors.append(_describe_error("data", e))
                raise
        return reader

    @classmethod
    def read_headers(cls, fo: str | os.PathLike | BinaryIO) -> "Reader":
        """Read only the headers -- no data -- of a biopac file."""
        with open_or_yield(fo, "rb") as io:
            reader = cls(io)
            try:
                reader._read_headers()
            except READ_EXCEPTIONS as e:
                # Consume the exception so the caller still gets the reader.
                reader._record_read_error("headers", e)
        return reader

    def _record_read_error(self, what, e):
        """Note an error we're recovering from, so it isn't silently lost."""
        message = _describe_error(what, e)
        logger.warning(message)
        self.read_errors.append(message)

    def stream(self, channel_indexes=None, target_chunk_size=CHUNK_SIZE):
        """Set up and retun an iterator for streaming data."""
        if self.datafile is None:
            self._read_headers()
        if self.is_compressed:
            raise TypeError("Streaming is not supported for compressed files")

        self.data_reader = DataReader(
            self.acq_file, self.datafile, self.data_start_offset
        )
        return self.data_reader.stream(channel_indexes, target_chunk_size)

    @property
    def is_compressed(self) -> bool:
        if self.graph_header is None:
            raise RuntimeError("Headers have not been read yet")
        return self.graph_header.compressed

    def __repr__(self):
        return f"Reader('{self.acq_file}')"

    # Header reading methods
    def _read_headers(self):
        # Read graph header and extract the necessary bits
        self.graph_header = GraphHeaderReader.bootstrap(self.acq_file)
        self.headers.append(self.graph_header)
        self.file_revision = self.graph_header.file_revision
        self.version_string = version_string_guess(self.file_revision)
        self.byte_order_char = self.graph_header.byte_order_char
        self.encoding = self.graph_header.encoding
        channel_count = self.graph_header.channel_count

        # We'll use this to read the rest of the headers
        self.header_reader = HeaderReader(
            self.acq_file, self.byte_order_char, self.file_revision, self.encoding
        )

        # Read padding headers
        pad_offset = self.graph_header.effective_len_bytes
        pad_headers = self.header_reader.multi_headers(
            self.graph_header.expected_padding_headers,
            pad_offset,
            bh.UnknownPaddingHeader,
        )
        self.headers.extend(pad_headers)

        # Read channel headers
        channel_offset = pad_offset + sum(
            [ph.effective_len_bytes for ph in pad_headers]
        )
        # skip past the unknown padding header

        channel_headers = self.header_reader.multi_headers(
            channel_count, channel_offset, bh.ChannelHeader
        )
        self.channel_headers = channel_headers
        self.headers.extend(channel_headers)

        for i, ch in enumerate(channel_headers):
            if logger.isEnabledFor(logging.DEBUG):
                logger.debug("Channel header %s: %s", i, ch.data)

        # Read foreign header. multi_headers() walked the channel headers by
        # each one's own effective_len_bytes, so the end of the last one is
        # where the foreign header starts -- don't re-derive it by assuming
        # every channel header is the same length.
        last_channel_header = channel_headers[-1]
        foreign_offset = (
            last_channel_header.offset + last_channel_header.effective_len_bytes
        )
        foreign_header = self.header_reader.single_header(
            foreign_offset, bh.ForeignHeader
        )
        self.headers.append(foreign_header)

        # Read channel dtype headers
        dtype_offset = foreign_offset + foreign_header.effective_len_bytes
        cdr = ChannelDTypeHeaderReader(self.header_reader)
        channel_dtype_headers = cdr.scan_for_dtype_headers(dtype_offset, channel_count)
        if len(channel_dtype_headers) == 0:
            raise ValueError("Can't find valid channel data type headers")
        self.channel_dtype_headers = channel_dtype_headers
        self.headers.extend(channel_dtype_headers)

        # We just read the data type headers, so we should be at the start of
        # the data
        self.data_start_offset = self.header_reader.acq_file.tell()

        logger.debug("Computed data start offset: %s", self.data_start_offset)

        # Calculate samples per second
        self.samples_per_second = 1000 / self.graph_header.sample_time

        # Create datafile
        logger.debug("Allocating a Datafile")
        self.datafile = Datafile(
            graph_header=self.graph_header,
            channel_headers=channel_headers,
            channel_dtype_headers=channel_dtype_headers,
            samples_per_second=self.samples_per_second,
        )
        if hasattr(self.acq_file, "name"):
            self.datafile.name = self.acq_file.name
        self.datafile.version_string = self.version_string

        self.data_length = self._data_length_bytes(
            channel_headers, channel_dtype_headers
        )
        logger.debug("Computed data length: %s", self.data_length)

        # In compressed files, markers come before compressed data. But
        # data_length is 0 for compressed files.
        marker_start_offset = self.data_start_offset + self.data_length
        self.marker_reader = MarkerReader.create_marker_reader(self.header_reader)
        self.event_markers = self.marker_reader.read_markers(
            marker_start_offset, self.graph_header.sample_time
        )
        self.headers.extend(self.marker_reader.all_headers)
        self.datafile.event_markers = self.event_markers

        # We should be right at the start of the journal data (if it exists)
        journal_offset = self.acq_file.tell()
        self.journal_reader = JournalReader.create_journal_reader(self.header_reader)
        try:
            self.journal = self.journal_reader.read_journal(journal_offset)
            self.datafile.journal = self.journal
            self.headers.extend(self.journal_reader.all_headers)
        except READ_EXCEPTIONS as e:
            self._record_read_error("journal", e)

        # Read compression headers if needed
        self._read_compression_headers_if_compressed(channel_count)

    def _data_length_bytes(self, channel_headers, channel_dtype_headers):
        if self.is_compressed:
            return 0
        return sum(
            [
                (ch.point_count * cdh.sample_size)
                # Lengths can differ on damaged files; read what we can
                for ch, cdh in zip(channel_headers, channel_dtype_headers, strict=False)
            ]
        )

    def _read_compression_headers_if_compressed(self, channel_count):
        # We need to have read the markers and journal; this puts us
        # at the correct file offset.
        if not self.is_compressed:
            return

        main_ch_offset = self.acq_file.tell()
        main_compression_header = self.header_reader.single_header(
            main_ch_offset, bh.MainCompressionHeader
        )
        self.headers.append(main_compression_header)
        cch_offset = main_ch_offset + main_compression_header.effective_len_bytes

        # read_data will need the compression headers
        self.channel_compression_headers = self.header_reader.multi_headers(
            channel_count, cch_offset, bh.ChannelCompressionHeader
        )
        if self.datafile is not None:
            self.datafile.channel_compression_headers = self.channel_compression_headers
        self.headers.extend(self.channel_compression_headers)

    def _read_data(self, channel_indexes, target_chunk_size=CHUNK_SIZE):
        self.data_reader = DataReader(
            self.acq_file, self.datafile, self.data_start_offset
        )

        self.data_reader.read_data(channel_indexes, target_chunk_size)


@contextmanager
def open_or_yield(thing, mode):
    """If 'thing' is a string or Path, open it and yield it. Otherwise, yield it.

    This lets you use a filename, open file, other IO object. If 'thing' was
    a filename, the file is guaranteed to be closed after yielding.
    """
    if isinstance(thing, IOBase):
        yield thing
    else:
        logger.debug("Opening file: %s", thing)
        with open(thing, mode) as f:
            yield (f)
