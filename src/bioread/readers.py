# coding: utf8
# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2010 Board of Regents of the University of Wisconsin System
#
# Written by John Ollinger <ollinger@wisc.edu> and Nate Vack <njvack@wisc.edu>
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison

from __future__ import with_statement
import struct

import numpy as np

from file_versions import *
from headers import GraphHeader, ChannelHeader, ChannelDTypeHeader
from headers import ForeignHeader
from biopac import Datafile, Channel

class AcqReader(object):
    """ 
    Main class for reading AcqKnowledge files. You'll probably call it like:
    >>> data = AcqReader.read("some_file.acq")
    """
    
    def __init__(self, acq_file):
        self.acq_file = acq_file
        # This must be set by _set_order_and_version
        self.byte_order_flag = None
        self.file_revision = None

    @classmethod
    def read(cls, filename):
        """
        The main access for this 
        """
        df = None
        with open(filename, 'rb') as f:
            reader = cls(f)
            return reader.read_file()
            
    def read_file(self):
        self.__setup()
        samples_per_second = 1000/self.graph_header['dSampleTime']
        df = Datafile(
            graph_header=self.graph_header, 
            channel_headers=self.channel_headers,
            foreign_header=self.foreign_header,
            channel_dtype_headers=self.channel_dtype_headers,
            samples_per_second=samples_per_second
            )
        self.channels = self.__build_channels(native_sps=samples_per_second)
        self.__read_data(self.channels)
        df.channels = self.channels
        self.data_file = df
        return self
        
    def __setup(self):
        if self.byte_order_flag is not None:
            return
        # TODO: Extract this into a factory class
        self.__set_order_and_version()
        self.__read_headers()
    
    def __read_headers(self):
        # Shorthand
        v = self.file_revision
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
    
    
    def __build_channels(self, native_sps=0.0):
        # Build empty channels, ready to get data from the file.
        channels = []
        # For building raw data arrays
        np_map = {
            1: np.float64,
            2: np.int16
        }
        fmt_map = {
            1: 'd',
            2: 'h'
        }
        
        for i in range(len(self.channel_headers)):
            ch = self.channel_headers[i].data
            cdh = self.channel_dtype_headers[i].data
            data = np.zeros(ch['lBufLength'], np_map[cdh['nType']])
            divider = ch['nVarSampleDivider']
            samples_per_second=float(native_sps)/divider
            chan = Channel(
                freq_divider=divider, raw_scale_factor=ch['dAmplScale'],
                raw_offset=ch['dAmplOffset'], raw_data=data,
                name=ch['szCommentText'], units=ch['szUnitsText'],
                fmt_str=fmt_map[cdh['nType']],
                samples_per_second=samples_per_second
            )
            channels.append(chan)
        return channels
    
    def __read_data(self, channels):
        # The data in the file are interleaved, so we'll potentially have
        # a different amount of data to read at each time slice.
        # It's possible we won't have any data for some time slices, I think.
        # The BIOPAC engineers tell you not to even try reading interleaved
        # data. Wusses.
        
        # This seems to be the same for all channels, but it's not specced.
        # This method should prevent us from leaving data from some channels.
        n_guesses = [c.freq_divider*c.raw_data.shape[0] for c in channels]
        max_n = max(n_guesses)
        
        self.acq_file.seek(self.data_start_offset)
        for i in xrange(max_n):
            sample_channels = [c for c in channels if i % c.freq_divider == 0]
            slice_fmt = self.byte_order_flag+''.join(
                [c.fmt_str for c in sample_channels])
            data = self.acq_file.read(struct.calcsize(slice_fmt))
            samples = struct.unpack(slice_fmt, data)
            for chan, samp in zip(sample_channels, samples):
                d_index = i//chan.freq_divider
                chan.raw_data[d_index] = samp
                
    
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
        self.file_revision = bp[0]
        
def pch(r):
    return [h.data for h in r.channel_dtype_headers]