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
import bioread.file_revisions as rev
from bioread import headers as bh
from bioread.biopac import EventMarker

logger = logging.getLogger("bioread")


class MarkerReader:
    """
    Base class for reading marker data from a BIOPAC file.
    """
    def __init__(self, header_reader):
        """
        Initialize a MarkerReader.
        
        Parameters
        ----------
        header_reader : HeaderReader
            The header reader to use for reading headers
        """
        self.header_reader = header_reader
        self.acq_file = header_reader.acq_file
        self.encoding = header_reader.encoding
        self.file_revision = header_reader.file_revision
        
    @staticmethod
    def create_marker_reader(header_reader, file_revision):
        """
        Factory method to create the appropriate marker reader based on file revision.
        
        Parameters
        ----------
        header_reader : HeaderReader
            The header reader to use
        file_revision : int
            The file revision
            
        Returns
        -------
        MarkerReader
            The appropriate marker reader
        """
        if file_revision >= rev.V_400B:
            return V4MarkerReader(header_reader)
        else:
            return V2MarkerReader(header_reader)
            
    def read_markers(self, marker_start_offset, graph_header):
        """
        Read the marker data.
        
        Parameters
        ----------
        marker_start_offset : int
            The offset where markers start
        graph_header : header
            The graph header
            
        Returns
        -------
        tuple
            (marker_header, marker_item_headers, event_markers)
        """
        raise NotImplementedError("Subclasses must implement read_markers")


class V2MarkerReader(MarkerReader):
    """
    Reader for version 2 marker data.
    """
    def read_markers(self, marker_start_offset, graph_header):
        """
        Read version 2 markers.
        
        Parameters
        ----------
        marker_start_offset : int
            The offset where markers start
        graph_header : header
            The graph header
            
        Returns
        -------
        tuple
            (marker_header, marker_item_headers, event_markers)
        """
        logger.debug("Reading markers starting at %s" % marker_start_offset)
        
        # Read marker header
        marker_header = self.header_reader.single_header(
            marker_start_offset, bh.V2MarkerHeader)
            
        # Read marker items
        marker_item_headers, event_markers = self._read_marker_items(
            marker_header, bh.V2MarkerItemHeader, graph_header)
            
        # Read marker metadata if needed
        marker_metadata_pre_header = None
        marker_metadata_headers = None
        if self.file_revision >= rev.V_381 and self.file_revision <= rev.V_400B:
            marker_metadata_pre_header, marker_metadata_headers = self._read_v2_marker_metadata(event_markers)
            
        return marker_header, marker_item_headers, event_markers, marker_metadata_pre_header, marker_metadata_headers
        
    def _read_marker_items(self, marker_header, marker_item_header_class, graph_header):
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
        
        for i in range(marker_header.marker_count):
            mih = self.header_reader.single_header(
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
        
    def _read_v2_marker_metadata(self, event_markers):
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
        marker_metadata_pre_header = self.header_reader.single_header(
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

        marker_metadata_headers = self.header_reader.multi_headers(
            marker_metadata_pre_header.item_count,
            self.acq_file.tell(),
            bh.V2MarkerMetadataHeader
        )
        for mh in marker_metadata_headers:
            event_markers[mh.marker_index].color = mh.rgba_color
            event_markers[mh.marker_index].tag = mh.marker_tag
            
        return marker_metadata_pre_header, marker_metadata_headers


class V4MarkerReader(MarkerReader):
    """
    Reader for version 4 marker data.
    """
    def read_markers(self, marker_start_offset, graph_header):
        """
        Read version 4 markers.
        
        Parameters
        ----------
        marker_start_offset : int
            The offset where markers start
        graph_header : header
            The graph header
            
        Returns
        -------
        tuple
            (marker_header, marker_item_headers, event_markers)
        """
        logger.debug("Reading markers starting at %s" % marker_start_offset)
        
        # Read marker header
        marker_header = self.header_reader.single_header(
            marker_start_offset, bh.V4MarkerHeader)
            
        # Read marker items
        marker_item_headers, event_markers = self._read_marker_items(
            marker_header, bh.V4MarkerItemHeader, graph_header)
            
        return marker_header, marker_item_headers, event_markers, None, None
        
    def _read_marker_items(self, marker_header, marker_item_header_class, graph_header):
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
        
        for i in range(marker_header.marker_count):
            mih = self.header_reader.single_header(
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