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
import struct

import bioread.file_revisions as rev
from bioread import headers as bh

# How far past the foreign data header we're willing to go looking for the
# channel dtype headers
MAX_DTYPE_SCANS = 4096

logger = logging.getLogger("bioread")


class HeaderReader:
    """
    A class to handle reading generic headers from a BIOPAC file.
    """
    def __init__(self, acq_file, byte_order_char=None, file_revision=None, encoding=None):
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
        
    def set_order_and_version(self):
        """
        Determine the byte order and file revision of the file.
        
        Returns
        -------
        tuple
            (file_revision, byte_order_char, encoding)
        """
        # Try unpacking the version string in both a bid and little-endian
        # fashion. Version string should be a small, positive integer.
        # It doesn't matter which graph header class we use; we're only using
        # the file revision field, which is the same for all graph headers
        graph_reads = [
            bh.GraphHeaderPre4(rev.V_ALL, bom)
            for bom in ['<', '>']
        ]
        for graph_header in graph_reads:
            graph_header.unpack_from_file(self.acq_file, 0)
            logger.debug(f"Interpreting file revision with byte order {graph_header.byte_order_char}: {graph_header.file_revision}")
            logger.debug(f"Graph header: {graph_header.raw_data.hex()}")
        
        rev_bom = [
            (graph_header._struct.lVersion, graph_header.byte_order_char)
            for graph_header in graph_reads
        ]
        rev_bom.sort()
        file_revision = rev_bom[0][0]
        byte_order_char = rev_bom[0][1]
        
        # Guess at file encoding -- I think that everything before acq4 is
        # in latin1 and everything newer is utf-8
        logger.debug(f"File revision: {file_revision}")
        logger.debug(f"Byte order: {byte_order_char}")
        if file_revision < rev.V_400B:
            encoding = 'latin1'
        else:
            encoding = 'utf-8'
            
        self.file_revision = file_revision
        self.byte_order_char = byte_order_char
        self.encoding = encoding
        
        return (file_revision, byte_order_char, encoding)
    
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
            logger.debug(f"Reading {h_class} at offset {h_offset}")
            h = h_class(self.file_revision,
                        self.byte_order_char,
                        encoding=self.encoding)
            try:
                h.unpack_from_file(self.acq_file, h_offset)
            except Exception as e:
                # If something goes wrong with reading, we want to log the
                # error and reraise it -- read_headers() and read_data() will
                # handle the error there.
                logger.error(
                    f"Error reading {h_class} at offset {h_offset}: {e}")
                raise e
            logger.debug(f"Read {h.struct_length} bytes: {h.data}")
            last_h_len = h.effective_len_bytes
            headers.append(h)
        return headers
    
    def scan_for_dtype_headers(self, start_index, channel_count):
        """
        Scan for channel dtype headers.
        
        Sometimes the channel dtype headers don't seem to be right after the
        foreign data header, and I can't find anything that directs me to the
        proper location.
        As a gross hack, we can scan forward until we find something
        potentially valid.
        
        Parameters
        ----------
        start_index : int
            The index to start scanning from
        channel_count : int
            The number of channels
            
        Returns
        -------
        list
            A list of channel dtype headers, or None if none were found
        data_start_offset : int
            The offset where the data starts
        """
        logger.debug('Scanning for start of channel dtype headers')
        for i in range(MAX_DTYPE_SCANS):
            dtype_headers = self.multi_headers(
                channel_count, start_index + i, bh.ChannelDTypeHeader)
            if all([h.possibly_valid for h in dtype_headers]):
                logger.debug(f"Found at {start_index + i}")
                data_start_offset = self.acq_file.tell()
                return dtype_headers, data_start_offset
        logger.warning(
            f"Couldn't find valid dtype headers, tried {MAX_DTYPE_SCANS} times"
        )
        return None, None
    
    def read_journal_v2(self):
        """
        Read a version 2 journal.
        
        Returns
        -------
        tuple
            (journal_header, journal)
        """
        logger.debug("Reading journal starting at %s" % self.acq_file.tell())
        logger.debug(self.acq_file.tell())
        journal_header = self.single_header(
            self.acq_file.tell(), bh.V2JournalHeader)
        if not journal_header.tag_value_matches_expected():
            raise ValueError(
                f"Journal header tag is {journal_header.tag_value_hex}, expected {bh.V2JournalHeader.EXPECTED_TAG_VALUE_HEX}"
            )
        journal = self.acq_file.read(
            journal_header.journal_len).decode(
                self.encoding, errors='ignore').strip('\0')
        return journal_header, journal
    
    def read_journal_v4(self):
        """
        Read a version 4 journal.
        
        Returns
        -------
        tuple
            (journal_header, journal, journal_length_header)
        """
        journal_length_header = self.single_header(
            self.acq_file.tell(),
            bh.V4JournalLengthHeader)
        journal_len = journal_length_header.journal_len
        journal = None
        journal_header = None
        jh = bh.V4JournalHeader(
            self.file_revision, self.byte_order_char)
        # If journal_length_header.journal_len is small, we don't have a
        # journal to read.
        if (jh.effective_len_bytes <= journal_len):
            journal_header = self.single_header(
                self.acq_file.tell(),
                bh.V4JournalHeader)
            logger.debug("Reading {0} bytes of journal at {1}".format(
                journal_header.journal_len,
                self.acq_file.tell()))
            journal = self.acq_file.read(
                journal_header.journal_len).decode(
                    self.encoding, errors='ignore').strip('\0')
        # Either way, we should seek to this point.
        self.acq_file.seek(journal_length_header.data_end)
        return journal_header, journal, journal_length_header
    
    def read_marker_items(self, marker_header, marker_item_header_class, graph_header):
        """
        Read marker items.
        
        Parameters
        ----------
        marker_header : header
            The marker header
        marker_item_header_class : class
            The marker item header class
        graph_header : header
            The graph header
            
        Returns
        -------
        tuple
            (marker_item_headers, event_markers)
        """
        event_markers = []
        marker_item_headers = []
        from bioread.biopac import EventMarker
        
        for i in range(marker_header.marker_count):
            mih = self.single_header(
                self.acq_file.tell(), marker_item_header_class)
            marker_text_bytes = self.acq_file.read(mih.text_length)
            marker_text = marker_text_bytes.decode(
                self.encoding, errors='ignore').strip('\0')
            marker_item_headers.append(mih)
            # We don't have the channel_order_map yet, so we'll just store None for now
            event_markers.append(EventMarker(
                time_index=(mih.sample_index * graph_header.sample_time) / 1000,
                sample_index=mih.sample_index,
                text=marker_text,
                channel_number=mih.channel_number,
                channel=None,
                date_created_ms=mih.date_created_ms,
                type_code=mih.type_code))
        
        return marker_item_headers, event_markers
    
    def read_v2_marker_metadata(self, event_markers):
        """
        Read version 2 marker metadata.
        
        Parameters
        ----------
        event_markers : list
            The event markers to update with metadata
            
        Returns
        -------
        tuple
            (marker_metadata_pre_header, marker_metadata_headers)
        """
        marker_metadata_pre_header = self.single_header(
            self.acq_file.tell(),
            bh.V2MarkerMetadataPreHeader
        )
        # Sometimes the marker metadata headers are not present -- if we see
        # that the first four bytes of the marker metadata preheader are
        # actually the start of the journal header (0x44332211), rewind the
        # file to the start of the marker metadata preheader and return
        if marker_metadata_pre_header.tag_value == bh.V2JournalHeader.EXPECTED_TAG_VALUE: 
            self.acq_file.seek(marker_metadata_pre_header.offset)
            logger.debug("No marker metadata headers found")
            return marker_metadata_pre_header, None

        marker_metadata_headers = self.multi_headers(
            marker_metadata_pre_header.item_count,
            self.acq_file.tell(),
            bh.V2MarkerMetadataHeader
        )
        for mh in marker_metadata_headers:
            event_markers[mh.marker_index].color = mh.rgba_color
            event_markers[mh.marker_index].tag = mh.marker_tag
            
        return marker_metadata_pre_header, marker_metadata_headers 