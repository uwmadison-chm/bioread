# coding: utf8
# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2010 Board of Regents of the University of Wisconsin System
#
# Written by John Ollinger <ollinger@wisc.edu> and Nate Vack <njvack@wisc.edu>
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison
# Project home: http://github.com/njvack/bioread

from __future__ import with_statement
import struct
import zlib

import numpy as np

from file_revisions import *
from headers import GraphHeader, ChannelHeader, ChannelDTypeHeader
from headers import ForeignHeader, MainCompressionHeader
from headers import ChannelCompressionHeader
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
        self.samples_per_second = None
        self.graph_header = None
        self.channel_headers = []
        self.foreign_header = None
        self.channel_dtype_headers = []
        self.main_compression_header = None
        self.channel_compression_headers = []
        self.data_start_offset = None

    @classmethod
    def read_file(cls, fo):
        """
        The main method to quickly read a biopac file into memory.

        fo: The name of the file to read, or a file-like object

        returns: biopac.Datafile
        """
        df = None
        if type(fo) == str:
            with open(fo, 'rb') as f:
                reader = cls(f)
                return reader.read()
        else:
            reader = cls(fo)
            return reader.read()

    def read(self):
        self._read_headers()
        self.samples_per_second = 1000/self.graph_header.sample_time
        df = Datafile(
            graph_header=self.graph_header,
            channel_headers=self.channel_headers,
            foreign_header=self.foreign_header,
            channel_dtype_headers=self.channel_dtype_headers,
            samples_per_second=self.samples_per_second)

        self._read_data()
        df.channels = self.channels
        self.data_file = df
        return self.data_file

    def _read_headers(self):
        if self.byte_order_flag is None:
            self.__set_order_and_version()
            
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
        if self.graph_header.compressed:
            self.__read_compression_headers()

    def __read_compression_headers(self):
        main_ch_start = self.data_start_offset
        self.main_compression_header = self.__single_header(main_ch_start,
            MainCompressionHeader)

        cch_start = (main_ch_start +
            self.main_compression_header.effective_len_bytes)
        self.channel_compression_headers = self.__multi_headers(
            self.graph_header.channel_count, cch_start,
            ChannelCompressionHeader)

    def __single_header(self, start_offset, h_class):
        return self.__multi_headers(1, start_offset, h_class)[0]

    def __multi_headers(self, num, start_offset, h_class):
        headers = []
        last_h_len = 0 # This will be changed when reading the channel headers
        h_offset = start_offset
        for i in range(num):
            h_offset += last_h_len
            h = h_class(self.file_revision, self.byte_order_flag)
            h.unpack_from_file(self.acq_file, h_offset)
            last_h_len = h.effective_len_bytes
            headers.append(h)
        return headers

    def __build_channels(self):
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
            chan_samp_per_sec=float(self.samples_per_second)/divider
            chan = Channel(
                freq_divider=divider, raw_scale_factor=ch.raw_scale,
                raw_offset=ch.raw_offset, raw_data=data,
                name=ch.name, units=ch.units,
                fmt_str=fmt_map[cdh.type_code],
                samples_per_second=chan_samp_per_sec)
            channels.append(chan)
        return channels

    def _read_data(self):
        self.channels = self.__build_channels()
        if self.graph_header.compressed:
            self.__read_data_compressed(self.channels)
        else:
            self.__read_data_uncompressed(self.channels)

    def __read_data_compressed(self, channels):
        # At least in post-4.0 files, the compressed data isn't interleaved at
        # all. It's stored in uniform compressed blocks -- this probably
        # compresses far better than interleaved data.
        # Strangely, the compressed data seems to always be little-endian.
        for i in range(len(channels)):
            cch = self.channel_compression_headers[i]
            chan = channels[i]
            # Data seems to be little-endian regardless of the rest of the file
            fmt_str = '<'+chan.fmt_str
            self.acq_file.seek(cch.compressed_data_offset)
            comp_data = self.acq_file.read(cch.compressed_data_len)
            decomp_data = zlib.decompress(comp_data)
            # raw_data starts out as zeros. Yeah, this feels hacky to me, too.
            np.add(chan.raw_data, np.fromstring(decomp_data, fmt_str),
                chan.raw_data)

    def __read_data_uncompressed(self, channels):
        # The data in the file are interleaved, so we'll potentially have
        # a different amount of data to read at each time slice.
        # It's possible we won't have any data for some time slices, I think.
        # The BIOPAC engineers tell you not to even try reading interleaved
        # data. Wusses.

        # This is, unfortunately, not necessarily the same for each channel.
        n_guesses = [c.freq_divider*c.raw_data.shape[0] for c in channels]
        max_n = max(n_guesses)
        self.acq_file.seek(self.data_start_offset)
        for i in xrange(max_n):
            sample_channels = [c for c in channels if c.should_sample_at(i)]
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
