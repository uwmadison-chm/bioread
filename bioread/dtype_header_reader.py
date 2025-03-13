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
from typing import Optional, Tuple, List
from bioread import headers as bh
from bioread.header_reader import HeaderReader

# How far past the foreign data header we're willing to go looking for the
# channel dtype headers
MAX_DTYPE_SCANS = 4096

logger = logging.getLogger("bioread")


class DTypeHeaderReader:
    """
    A class to handle reading channel data type headers from a BIOPAC file.
    """

    def __init__(self, header_reader: HeaderReader) -> None:
        """
        Initialize a DTypeHeaderReader.

        Parameters
        ----------
        header_reader : HeaderReader
            The header reader to use for reading headers
        """
        self.header_reader = header_reader
        self.acq_file = header_reader.acq_file

    def scan_for_dtype_headers(
        self, start_index: int, channel_count: int
    ) -> Tuple[Optional[List[bh.ChannelDTypeHeader]], Optional[int]]:
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
        tuple
            (dtype_headers, data_start_offset) or (None, None) if not found
        """
        logger.debug("Scanning for start of channel dtype headers")
        for i in range(MAX_DTYPE_SCANS):
            dtype_headers = self.header_reader.multi_headers(
                channel_count, start_index + i, bh.ChannelDTypeHeader
            )
            if all(h.possibly_valid for h in dtype_headers):
                logger.debug(f"Found at {start_index + i}")
                data_start_offset = self.acq_file.tell()
                return dtype_headers, data_start_offset
        logger.warning(
            f"Couldn't find valid dtype headers, tried {MAX_DTYPE_SCANS} times"
        )
        return None, None
