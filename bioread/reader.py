# coding: utf8
# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2025 Board of Regents of the University of Wisconsin System
#
# Written Nate Vack <njvack@wisc.edu> with research from John Ollinger
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison
# Extended by Alexander Schlemmer.

import struct
from contextlib import contextmanager
from io import IOBase

import numpy as np

import bioread.file_revisions as rev
from bioread import headers as bh
from bioread.biopac import Datafile
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


class Reader:
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
            header_and_journal = self.journal_reader.read_journal()
            if header_and_journal is None:
                return
            self.journal_header, self.journal = header_and_journal
        except READ_EXCEPTIONS as e:
            logger.error(f"Error reading journal: {e}")
            self.read_errors.append(f"Error reading journal: {e}")
            raise
            
        self.datafile.journal_header = self.journal_header
        self.datafile.journal = self.journal

    # Marker reading methods
    def _read_markers(self):
        
        self.marker_reader.read_markers(self.marker_start_offset, self.graph_header.sample_time)
        self.event_markers = self.marker_reader.event_markers
        self.marker_header = self.marker_reader.marker_header
        self.marker_item_headers = self.marker_reader.marker_item_headers
        self.marker_metadata_pre_header = self.marker_reader.marker_metadata_pre_header
        self.marker_metadata_headers = self.marker_reader.marker_metadata_headers
            
        # Update channel references in event markers
        for marker in self.event_markers:
            marker.channel = self.datafile.channel_order_map.get(
                marker.channel_number)
                
        self.datafile.event_markers = self.event_markers
        self.datafile.marker_header = self.marker_header
        self.datafile.marker_item_headers = self.marker_item_headers

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
    """ If 'thing' is a string or Path, open it and yield it. Otherwise, yield it.

    This lets you use a filename, open file, other IO object. If 'thing' was
    a filename, the file is guaranteed to be closed after yielding.
    """
    if isinstance(thing, IOBase):
        yield thing
    else:
        logger.debug(f"Opening file: {thing}")
        with open(thing, mode) as f:
            yield(f)
