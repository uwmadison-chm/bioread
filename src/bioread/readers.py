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

from file_revisions import *
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
    def read_file(cls, filename):
        """
        The main method to quickly read a biopac file into memory.

        filename: The name of the file to read.

        returns: biopac.Datafile
        """
        df = None
        with open(filename, 'rb') as f:
            reader = cls(f)
            return reader.read()

    def read(self):
        self.__setup()
        samples_per_second = 1000/self.graph_header.sample_time
        df = Datafile(
            graph_header=self.graph_header,
            channel_headers=self.channel_headers,
            foreign_header=self.foreign_header,
            channel_dtype_headers=self.channel_dtype_headers,
            samples_per_second=samples_per_second)

        self.channels = self.__build_channels(native_sps=samples_per_second)
        self.__read_data(self.channels)
        df.channels = self.channels
        self.data_file = df
        return self.data_file

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
        self.graph_header = self.__single_header(0, GraphHeader)
        channel_count = self.graph_header.channel_count

        ch_start = self.graph_header.effective_len_bytes
        self.channel_headers = self.__multi_headers(channel_count, 
            ch_start, ChannelHeader)
        ch_len = self.channel_headers[0].effective_len_bytes

        fh_start = ch_start + len(self.channel_headers)*ch_len
        self.foreign_header = self.__single_header(fh_start, ForeignHeader)

        cdh_start = fh_start + self.foreign_header.effective_len_bytes
        self.channel_dtype_headers = self.__multi_headers(channel_count, 
            cdh_start, ChannelDTypeHeader)
        cdh_len = self.channel_dtype_headers[0].effective_len_bytes
        
        self.data_start_offset = (cdh_start + (cdh_len * channel_count))
    
    def __single_header(self, start_offset, h_class):
        h = h_class(self.file_revision, self.byte_order_flag)
        h.unpack_from_file(self.acq_file, start_offset)
        return h
    
    def __multi_headers(self, num, start_offset, h_class):
        headers = []
        h_len = 0 # This will be changed when reading the channel headers
        for i in range(num):
            # OK for ch_len to be 0 on first iter
            h_offset = start_offset + i*h_len 
            h = h_class(self.file_revision, self.byte_order_flag)
            h.unpack_from_file(self.acq_file, h_offset)
            h_len = h.effective_len_bytes
            headers.append(h)
        return headers

    def __build_channels(self, native_sps=0.0):
        # Build empty channels, ready to get data from the file.
        channels = []
        # For building raw data arrays
        np_map = {
            1: np.float64,
            2: np.int16}
        fmt_map = {
            1: 'd',
            2: 'h'}

        for i in range(len(self.channel_headers)):
            ch = self.channel_headers[i]
            cdh = self.channel_dtype_headers[i]
            data = np.zeros(ch.point_count, np_map[cdh.type_code])
            divider = ch.frequency_divider
            samples_per_second=float(native_sps)/divider
            chan = Channel(
                freq_divider=divider, raw_scale_factor=ch.raw_scale,
                raw_offset=ch.raw_offset, raw_data=data,
                name=ch.name, units=ch.units,
                fmt_str=fmt_map[cdh.type_code],
                samples_per_second=samples_per_second)
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
                bof in byte_order_flags]

        # Limit to positive numbers, choose smallest.
        byte_order_versions = sorted([
            bp for bp in byte_order_versions if bp[0] > 0])
        bp = byte_order_versions[0]

        self.byte_order_flag = bp[1]
        self.file_revision = bp[0]


# Notes on compressed files:
# At data_start_offset, there's a long, containing the length of some header H1.
# At data_start_offset + len(H1), there's something.
# From there, it's 95 bytes to the start of the first channel header.
# Channel headers are 59+l_channel+l_unit bytes long.
# The start values' uses are unknown, but:
# bytes 43-46 are length of channel label (l_channel)
# bytes 47-50 are length of unit label (l_unit)
# bytes 51-54 are (related to) number of data points -- maybe size of 
#                 uncompressed data?
# bytes 55-58 are length of compressed data (l_comp)
# So, the compressed data for the first channel should start at offset:
# data_start_offset + 95+59+l_channel+l_unit
# and it should be l_comp bytes long.
# The next channel header follows immediately.
