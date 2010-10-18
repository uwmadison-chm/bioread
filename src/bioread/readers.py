# coding: utf8
# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2010 Board of Regents of the University of Wisconsin System
#
# Written by John Ollinger <ollinger@wisc.edu> and Nate Vack <njvack@wisc.edu>
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison

import struct

from struct_dict import StructDict
from headers import Header
from file_versions import *


class VersionedReader(object):
    def __init__(self, acq_file=None, byte_order_flag=None, file_version=None):
        self.acq_file = acq_file
        self.byte_order_flag = byte_order_flag
        self.file_version = file_version

    def read(self):
        """ 
        Right now, do the simplest thing that'll create the data structures
        we want. Soon, create a proper biopac.Datafile structure with Channels
        and such.
        """
        self.graph_header = self.__read_header(self.graph_header_sd, 0)
        self.channel_headers = []
        num_channels = self.graph_header['nChannels']
        graph_header_len = self.graph_header['lExtItemHeaderLen']
        channel_head_len = 0 # This will get changed!
        for i in range(num_channels):
            channel_offset = graph_header_len + i*channel_head_len
            channel_header = self.__read_header(
                self.channel_header_sd, channel_offset)
            # This will be the same for all channels; easier to set it for all
            channel_head_len = channel_header['lChanHeaderLen']
            self.channel_headers.append(channel_header)

        foreign_header_offset = (
            graph_header_len + num_channels*channel_head_len)
        self.foreign_header = self.__read_header(
            self.foreign_header_sd, foreign_header_offset)
        foreign_header_len = self.foreign_header['nLength']

        channel_dtype_offset = foreign_header_offset + foreign_header_len
        channel_dtype_len = self.channel_dtype_header_sd.len_bytes
        self.channel_dtype_headers = []
        for i in range(num_channels):
            offset = channel_dtype_offset + i*channel_dtype_len
            dt_header = self.__read_header(
                self.channel_dtype_header_sd, offset
            )
            self.channel_dtype_headers.append(dt_header)

        self.data_start_offset = (
            channel_dtype_offset + num_channels * channel_dtype_len)

    
    def __read_header(self, struct_dict, offset):
        h = Header(struct_dict)
        h.unpack_from_file(self.acq_file, offset)
        return h        

class AcqReader(VersionedReader):
    """ 
    Main class for reading AcqKnowledge files. You'll probably call it like:
    >>> data = AcqReader.read("some_file.acq")
    """
    
    def __init__(self, acq_file):
        self.acq_file = acq_file
        # This must be set by _set_order_and_version
        self.byte_order_flag = None

    def __setup(self):
        if self.byte_order_flag is not None:
            return
        # TODO: Extract this into a factory class
        self.__set_order_and_version()

    def __set_order_and_version(self):
        # Try unpacking the version string in both a bid and little-endian
        # fashion. Version string should be a small, positive integer.
        self.acq_file.seek(0)
        mh = self._graph_headers_for(V_ALL)
        ver_fmt_str = ''.join([s[1] for s in mh])
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