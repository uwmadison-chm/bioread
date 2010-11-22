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
from utils import lcm


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
                fmt_str=self.byte_order_flag+fmt_map[cdh.type_code],
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
            fmt_str = '<'+chan.fmt_str[1]
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

        # Using adapted algorithm by Sven Marnarch from:
        # http://stackoverflow.com/questions/4227990
        stream_byte_indexes = self.__stream_byte_indexes(channels)
        self.stream_byte_indexes = stream_byte_indexes
        self.block_len = len(stream_byte_indexes)
        # Allocate memory for our data stream -- it'll be padded if recording
        # stops in the middle of a block.
        self.data_len = sum([c.data_length for c in channels])
        self.num_blocks = int(np.ceil(float(self.data_len)/self.block_len))
        self.buf_len = self.block_len*self.num_blocks
        self.buf = np.zeros(self.buf_len, dtype=np.ubyte)
        self.buf[0:self.data_len] += np.fromfile(
            self.acq_file, np.ubyte, self.data_len)
        self.acq_file.seek(self.data_start_offset)
        # Now, partition the data into chunks of block_len
        self.buf = self.buf.reshape(-1, self.block_len)
        # and fill in the data.
        for i, ch in enumerate(channels):
            tmp = self.buf[:,stream_byte_indexes == i].ravel()
            tmp.dtype = ch.fmt_str
            np.add(tmp[0:ch.point_count], 0, ch.raw_data)

    def __stream_sample_indexes(self, channels):
        """
        Returns the shortest repeating pattern of samples that'll appear in
        our data stream. If our freq_dividers look like [1,2,4], we'll return
        [0,1,2,0,0,1,0]
        """
        dividers =[c.freq_divider for c in channels]
        channel_lcm = lcm(*dividers)
        # Make a list like [0,1,2,0,0,1,0]
        stream_sample_indexes = [
            ch_idx for pat_idx in range(channel_lcm)
            for ch_idx, div in enumerate(dividers)
            if pat_idx % div == 0]
        return stream_sample_indexes

    def __stream_byte_indexes(self, channels):
        """
        Returns the shortest repeating pattern of bytes that'll appear in
        our data stream. If our freq_dividers look like [1,2,4], and our
        sample_lengths are [2,2,8], we'll return:
        [0,0,1,1,2,2,2,2,2,2,2,2,0,0,0,0,1,1,0,0]
        """
        ssi = self.__stream_sample_indexes(channels)
        repeats = [channels[i].sample_size for i in ssi]
        return np.array(ssi).repeat(repeats)

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
