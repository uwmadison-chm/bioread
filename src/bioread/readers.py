# coding: utf8
# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2010 Board of Regents of the University of Wisconsin System
#
# Written by John Ollinger <ollinger@wisc.edu> and Nate Vack <njvack@wisc.edu>
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison

import struct

from file_versions import *
from headers import GraphHeader, ChannelHeader, ChannelDTypeHeader
from headers import ForeignHeader


class AcqReader(object):
    """ 
    Main class for reading AcqKnowledge files. You'll probably call it like:
    >>> data = AcqReader.read("some_file.acq")
    """
    
    def __init__(self, acq_file):
        self.acq_file = acq_file
        # This must be set by _set_order_and_version
        self.byte_order_flag = None
        self.file_version = None

    def read(self):
        self.__setup()
    
    def __setup(self):
        if self.byte_order_flag is not None:
            return
        # TODO: Extract this into a factory class
        self.__set_order_and_version()
        self.__read_headers()
    
    def __read_headers(self):
        # Shorthand
        v = self.file_version
        bof = self.byte_order_flag
        self.graph_header = GraphHeader(v, bof)
        self.graph_header.unpack_from_file(self.acq_file, 0)
        channel_count = self.graph_header['nChannels']
        
        gh_len = self.graph_header.effective_len_bytes
        ch_len = 0 # This will be changed when reading the channel headers
        self.channel_headers = []
        for i in range(channel_count):
            ch_offset = gh_len + i*ch_len # OK for ch_len to be 0 on first iter
            ch = ChannelHeader(v, bof)
            ch.unpack_from_file(self.acq_file, ch_offset)
            ch_len = ch.effective_len_bytes
            self.channel_headers.append(ch)
        
        fh_offset = gh_len + len(self.channel_headers)*ch_len
        self.foreign_header = ForeignHeader(v, bof)
        self.foreign_header.unpack_from_file(self.acq_file, fh_offset)
        
        cdh_len = 0 # Gets changed just like ch_len
        self.channel_dtype_headers = []
        for i in range(channel_count):
            cdh_offset = (fh_offset + self.foreign_header.effective_len_bytes +
                (i * cdh_len))
            cdh = ChannelDTypeHeader(v, bof)
            cdh.unpack_from_file(self.acq_file, cdh_offset)
            cdh_len = cdh.effective_len_bytes
            self.channel_dtype_headers.append(cdh)
        
        self.data_start_offset = (
            fh_offset + self.foreign_header.effective_len_bytes + 
            (cdh_len * channel_count)
        )
    
    def __set_order_and_version(self):
        # Try unpacking the version string in both a bid and little-endian
        # fashion. Version string should be a small, positive integer.
        self.acq_file.seek(0)
        # No byte order flag -- we're gonna figure it out.
        gh = GraphHeader(V_ALL, '') 
        ver_fmt_str = gh.format_string
        ver_len = struct.calcsize('<'+ver_fmt_str)
        ver_data = self.acq_file.read(ver_len)
        
        byte_order_flags = ['<', '>']
        # Try both ways.
        byte_order_versions = [
            (struct.unpack(bof+ver_fmt_str, ver_data)[1], bof) for 
                bof in byte_order_flags
        ]
        
        # Limit to positive numbers, choose smallest.
        byte_order_versions = sorted([
            bp for bp in byte_order_versions if bp[0] > 0
        ])
        bp = byte_order_versions[0]
        
        self.byte_order_flag = bp[1]
        self.file_version = bp[0]