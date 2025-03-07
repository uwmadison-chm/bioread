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
        self.encoding = header_reader.encoding
        
    @staticmethod
    def create_journal_reader(header_reader, file_revision):
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
        if file_revision <= rev.V_400B:
            return V2JournalReader(header_reader)
        else:
            return V4JournalReader(header_reader)
            
    def read_journal(self):
        """
        Read the journal data.
        
        Returns
        -------
        tuple
            (journal_header, journal)
        """
        raise NotImplementedError("Subclasses must implement read_journal")


class V2JournalReader(JournalReader):
    """
    Reader for version 2 journal data.
    """
    def read_journal(self):
        """
        Read a version 2 journal.
        
        Returns
        -------
        tuple
            (journal_header, journal)
        """
        logger.debug("Reading journal starting at %s" % self.acq_file.tell())
        logger.debug(self.acq_file.tell())
        journal_header = self.header_reader.single_header(
            self.acq_file.tell(), bh.V2JournalHeader)
        if not journal_header.tag_value_matches_expected():
            raise ValueError(
                f"Journal header tag is {journal_header.tag_value_hex}, expected {bh.V2JournalHeader.EXPECTED_TAG_VALUE_HEX}"
            )
        journal = self.acq_file.read(
            journal_header.journal_len).decode(
                self.encoding, errors='ignore').strip('\0')
        return journal_header, journal


class V4JournalReader(JournalReader):
    """
    Reader for version 4 journal data.
    """
    def read_journal(self):
        """
        Read a version 4 journal.
        
        Returns
        -------
        tuple
            (journal_header, journal, journal_length_header)
        """
        journal_length_header = self.header_reader.single_header(
            self.acq_file.tell(),
            bh.V4JournalLengthHeader)
        journal_len = journal_length_header.journal_len
        journal = None
        journal_header = None
        jh = bh.V4JournalHeader(
            self.header_reader.file_revision, self.header_reader.byte_order_char)
        # If journal_length_header.journal_len is small, we don't have a
        # journal to read.
        if (jh.effective_len_bytes <= journal_len):
            journal_header = self.header_reader.single_header(
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