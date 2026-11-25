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

        self.marker_header = None
        self.marker_item_headers = None
        self.event_markers = None
        self.marker_metadata_pre_header = None
        self.marker_metadata_headers = None
        self.all_headers = []

    @staticmethod
    def create_marker_reader(header_reader):
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
        if header_reader.file_revision >= rev.V_400B:
            return V4MarkerReader(header_reader)
        return V2MarkerReader(header_reader)

    def read_markers(self, marker_start_offset, sample_time):
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

        event_markers: List[EventMarker]
        """
        raise NotImplementedError("Subclasses must implement read_markers")


class V2MarkerReader(MarkerReader):
    """
    Reader for version 2 marker data.
    """

    def read_markers(self, offset, sample_time):
        """
        Read version 2 markers.

        Parameters
        ----------
        offset : int
            The offset where markers start

        sample_time : float

                Returns
        -------
        event_markers
        """
        logger.debug("Reading markers starting at %s" % offset)

        # Read marker header
        marker_header = self.header_reader.single_header(offset, bh.MarkerHeader)
        self.all_headers.append(marker_header)

        # Read marker items
        marker_item_offset = marker_header.offset + marker_header.effective_len_bytes
        event_markers = self._read_marker_items(
            marker_header.marker_count, marker_item_offset, sample_time
        )

        metadata_offset = self.acq_file.tell()

        if self.file_revision >= rev.V_381:
            # This consumes data and modifies event_markers and self.all_headers
            self._read_marker_metadata(event_markers, metadata_offset)

        return event_markers

    def _read_marker_items(self, marker_count, offset, sample_time):
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
        marker_item_headers = self.header_reader.multi_headers(
            marker_count, offset, bh.MarkerItemHeader
        )
        self.all_headers.extend(marker_item_headers)
        event_markers = []
        for mih in marker_item_headers:
            self.acq_file.seek(mih.offset + mih.struct_length)
            marker_text_bytes = self.acq_file.read(mih.text_length)
            marker_text = marker_text_bytes.decode(
                self.encoding, errors="ignore"
            ).strip("\0")
            logger.debug(f"Marker text: {marker_text}")
            event_markers.append(
                EventMarker(
                    time_index=(mih.sample_index * sample_time) / 1000,
                    sample_index=mih.sample_index,
                    text=marker_text,
                    channel_number=mih.channel_number,
                    channel=None,
                    date_created_ms=mih.date_created_ms,
                    type_code=mih.type_code,
                )
            )
        return event_markers

    def _read_marker_metadata(self, event_markers, metadata_offset):
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
            metadata_offset, bh.MarkerPreItemMetadataHeaderPre4
        )

        # Sometimes the marker metadata headers are not present -- if we see
        # that the first four bytes of the marker metadata preheader are
        # actually the start of the journal header (0x44332211), rewind the
        # file to the start of the marker metadata preheader and return
        if (
            marker_metadata_pre_header.tag_value
            == bh.JournalHeaderPre4.EXPECTED_TAG_VALUE
        ):
            self.acq_file.seek(marker_metadata_pre_header.offset)
            logger.debug("No marker metadata headers found")
            return
        self.all_headers.append(marker_metadata_pre_header)

        # Read marker metadata headers
        marker_metadata_headers = self.header_reader.multi_headers(
            len(event_markers), self.acq_file.tell(), bh.MarkerItemMetadataHeader
        )
        for mh in marker_metadata_headers:
            event_markers[mh.marker_index].color = mh.rgba_color
            event_markers[mh.marker_index].tag = mh.marker_tag
        self.all_headers.extend(marker_metadata_headers)


class V4MarkerReader(MarkerReader):
    """
    Reader for version 4 marker data.
    """

    def read_markers(self, offset, sample_time):
        """
        Read version 4 markers.

        Parameters
        ----------
        offset : int
            The offset where markers start
        sample_time : float
            The sample time

        Returns
        -------
        tuple
            (marker_header, marker_item_headers, event_markers)
        """
        logger.debug(f"Reading markers starting at {offset}")

        # Read marker header
        marker_header = self.header_reader.single_header(offset, bh.MarkerHeader)
        self.all_headers.append(marker_header)

        # Read marker items
        event_markers_offset = marker_header.offset + marker_header.effective_len_bytes
        event_markers = self._read_marker_items(
            marker_header.marker_count, event_markers_offset, sample_time
        )

        return event_markers

    def _read_marker_items(self, marker_count, offset, sample_time):
        """
        Read marker items.

        Parameters
        ----------
        marker_header : header
            The marker header
        marker_item_header_class : class
            The marker item header class
        sample_time : float
            The sample time

        Returns
        -------
        tuple
            (marker_item_headers, event_markers)
        """
        event_markers = []
        marker_item_headers = self.header_reader.multi_headers(
            marker_count, offset, bh.MarkerItemHeader
        )
        self.all_headers.extend(marker_item_headers)

        for mih in marker_item_headers:
            self.acq_file.seek(mih.offset + mih.struct_length)
            marker_text_bytes = self.acq_file.read(mih.text_length)
            marker_text = marker_text_bytes.decode(
                self.encoding, errors="ignore"
            ).strip("\0")
            logger.debug(f"Marker text: {marker_text}")
            # We don't have the channel_order_map yet, so we'll just store None for now
            event_markers.append(
                EventMarker(
                    time_index=(mih.sample_index * sample_time) / 1000,
                    sample_index=mih.sample_index,
                    text=marker_text,
                    channel_number=mih.channel_number,
                    channel=None,
                    date_created_ms=mih.date_created_ms,
                    type_code=mih.type_code,
                )
            )

        return event_markers
