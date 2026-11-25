# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2025 Board of Regents of the University of Wisconsin System
#
# Written Nate Vack <njvack@wisc.edu> with research from John Ollinger
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison
# Extended by Alexander Schlemmer.

import logging

import bioread.file_revisions as rev
from bioread import headers as bh

logger = logging.getLogger(__name__)


class HeaderReader:
    """
    A class to handle reading generic headers from a BIOPAC file.
    """

    def __init__(
        self, acq_file, byte_order_char=None, file_revision=None, encoding=None
    ):
        """
        Initialize a HeaderReader.

        Parameters
        ----------
        acq_file : file-like object
            The file to read headers from
        byte_order_char : str, optional
            The byte order character ('<' or '>')
        file_revision : int, optional
            The file revision
        encoding : str, optional
            The text encoding to use
        """
        self.acq_file = acq_file
        self.byte_order_char = byte_order_char
        self.file_revision = file_revision
        self.encoding = encoding

    def single_header(self, start_offset, h_class):
        """
        Read a single header from the file.

        Parameters
        ----------
        start_offset : int
            The offset to start reading from
        h_class : class
            The header class to use

        Returns
        -------
        header
            The header object
        """
        return self.multi_headers(1, start_offset, h_class)[0]

    def multi_headers(self, num, start_offset, h_class):
        """
        Read multiple headers from the file.

        Parameters
        ----------
        num : int
            The number of headers to read
        start_offset : int
            The offset to start reading from
        h_class : class
            The header class to use

        Returns
        -------
        list
            A list of header objects
        """
        headers = []
        last_h_len = 0  # This will be changed reading the channel headers
        h_offset = start_offset
        for i in range(num):
            h_offset += last_h_len
            logger.debug("Reading %s at offset %s", h_class.__name__, h_offset)
            h = h_class.for_revision(
                self.file_revision, self.byte_order_char, encoding=self.encoding
            )
            try:
                h.unpack_from_file(self.acq_file, h_offset)
            except Exception as e:
                # Say where it went wrong; Reader decides whether to recover.
                e.add_note(f"While reading {h_class.__name__} at offset {h_offset}")
                raise
            if logger.isEnabledFor(logging.DEBUG):
                # h.data builds a dict of every field, so only do it if needed
                logger.debug("Read %s bytes: %s", h.struct_length, h.data)
            last_h_len = h.effective_len_bytes
            headers.append(h)
        return headers


class GraphHeaderReader:
    @classmethod
    def bootstrap(cls, acq_file):
        """
        Determine the byte order and revision, and encoding of the file, and return
        a GraphHeader.
        """
        # Try unpacking the version string in both a bid and little-endian
        # fashion. Version string should be a small, positive integer.
        # It doesn't matter which graph header class we use; we're only using
        # the file revision field, which is the same for all graph headers
        graph_reads = [
            bh.GraphHeader.for_revision(rev.V_ALL, bom) for bom in ["<", ">"]
        ]
        for graph_header in graph_reads:
            graph_header.unpack_from_file(acq_file, 0)
            logger.debug(
                "Interpreting file revision with byte order %s: %s",
                graph_header.byte_order_char,
                graph_header.file_revision,
            )

        rev_bom = [
            (graph_header._struct.lVersion, graph_header.byte_order_char)
            for graph_header in graph_reads
        ]
        rev_bom.sort()
        file_revision = rev_bom[0][0]
        byte_order_char = rev_bom[0][1]

        # Guess at file encoding -- I think that everything before acq4 is
        # in latin1 and everything newer is utf-8
        logger.debug("File revision: %s", file_revision)
        logger.debug("Byte order: %s", byte_order_char)
        if file_revision < rev.V_400B:
            encoding = "latin1"
        else:
            encoding = "utf-8"

        graph_header = bh.GraphHeader.for_revision(
            file_revision, byte_order_char, encoding
        )
        graph_header.unpack_from_file(acq_file, 0)

        return graph_header


class ChannelDTypeHeaderReader:
    MAX_DTYPE_SCANS = 4096

    def __init__(self, header_reader):
        self.header_reader = header_reader

    def scan_for_dtype_headers(self, start_index, channel_count):
        """
        Scan for channel dtype headers.
        """
        logger.debug("Scanning for start of channel dtype headers")
        for i in range(self.MAX_DTYPE_SCANS):
            offset = start_index + i
            logger.debug("Scanning at offset %s", offset)
            dtype_headers = self.header_reader.multi_headers(
                channel_count, offset, bh.ChannelDTypeHeader
            )
            if all([h.possibly_valid for h in dtype_headers]):
                logger.debug("Found at %s", offset)
                return dtype_headers
        logger.warning(
            "Couldn't find valid dtype headers, tried %s times", self.MAX_DTYPE_SCANS
        )
        return []
