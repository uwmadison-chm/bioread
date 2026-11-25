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

logger = logging.getLogger("bioread")


class JournalReader:
    """
    Base class for reading journal data from a BIOPAC file.
    """

    def __init__(self, header_reader):
        """
        Initialize a JournalReader.

        Parameters
        ----------
        header_reader : HeaderReader
            The header reader to use for reading headers
        """
        self.header_reader = header_reader
        self.acq_file = header_reader.acq_file
        self.file_revision = header_reader.file_revision
        self.encoding = header_reader.encoding
        self.all_headers = []

    @staticmethod
    def create_journal_reader(header_reader):
        """
        Factory method to create the appropriate journal reader based on file revision.

        Parameters
        ----------
        header_reader : HeaderReader
            The header reader to use
        file_revision : int
            The file revision

        Returns
        -------
        JournalReader
            The appropriate journal reader
        """
        if header_reader.file_revision >= rev.V_400B:
            return V4JournalReader(header_reader)
        return V2JournalReader(header_reader)

    def read_journal(self, offset):
        """
        Read the journal data.

        Returns
        -------
        journal
        """
        raise NotImplementedError("Subclasses must implement read_journal")


class V2JournalReader(JournalReader):
    """
    Reader for version 2 journal data.
    """

    def read_journal(self, offset):
        """
        Read a version 2 journal.

        Returns
        -------
        journal: str
        """
        if self.header_reader.file_revision < rev.V_370:
            # We don't know how to read journals from before V_370
            logger.info(
                f"Can't read journals from file revision {self.header_reader.file_revision}"
            )
            return None
        logger.debug("Reading journal starting at %s" % self.acq_file.tell())
        journal_header = self.header_reader.single_header(offset, bh.JournalHeader)
        if not journal_header.tag_value_matches_expected():
            raise ValueError(
                f"Journal header tag is {journal_header.tag_value_hex}, expected {bh.JournalHeader.EXPECTED_TAG_VALUE_HEX}"
            )
        self.all_headers.append(journal_header)
        logger.debug(
            f"Reading {journal_header.journal_len} bytes of journal at {self.acq_file.tell()}"
        )
        journal = (
            self.acq_file.read(journal_header.journal_len)
            .decode(self.encoding, errors="ignore")
            .strip("\0")
        )
        return journal


class V4JournalReader(JournalReader):
    """
    Reader for version 4 journal data.
    """

    def read_journal(self, offset):
        """
        Read a version 4 journal.

        Returns
        -------
        journal: str
        """
        journal_length_header = self.header_reader.single_header(
            offset, bh.JournalLengthHeader
        )
        journal_length = journal_length_header.journal_len
        self.all_headers.append(journal_length_header)
        journal_offset = offset + journal_length_header.effective_len_bytes

        # We're just using this to get the journal length, it won't have data
        jh = bh.JournalHeader.for_revision(
            self.header_reader.file_revision, self.header_reader.byte_order_char
        )
        # If the journal length as reported by the journal length header is
        # less than the length of the journal header, we don't have a journal.
        journal = None
        if jh.effective_len_bytes <= journal_length:
            journal_header = self.header_reader.single_header(
                journal_offset, bh.JournalHeader
            )
            self.all_headers.append(journal_header)
            logger.debug(
                "Reading {0} bytes of journal at {1}".format(
                    journal_header.journal_len, self.acq_file.tell()
                )
            )
            journal = (
                self.acq_file.read(journal_header.journal_len)
                .decode(self.encoding, errors="ignore")
                .strip("\0")
            )

        # Either way, journal_length_header.data_end is the end of all the
        # journal data, so we seek to that point.
        self.acq_file.seek(journal_length_header.data_end)
        return journal
