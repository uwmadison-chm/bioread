# coding: utf8
# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2016 Board of Regents of the University of Wisconsin System
#
# Written Nate Vack <njvack@wisc.edu> with research from John Ollinger
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison
# Project home: http://github.com/njvack/bioread

from __future__ import division
import numpy as np


class Datafile(object):
    """
    A data file for the AcqKnowledge system. Generally, gonna be created
    from a file by reader.Reader.
    """

    def __init__(
            self,
            graph_header=None, channel_headers=None, foreign_header=None,
            channel_dtype_headers=None, samples_per_second=None, name=None,
            marker_header=None, marker_item_headers=None,):
        self.graph_header = graph_header
        self.channel_headers = channel_headers
        self.foreign_header = foreign_header
        self.channel_dtype_headers = channel_dtype_headers
        self.samples_per_second = samples_per_second
        self.name = name
        self.marker_header = marker_header
        self.marker_item_headers = marker_item_headers
        self.markers = None
        self.journal_header = None
        self.journal = None
        self.__named_channels = None
        self.__time_index = None

        self.channels = self.__build_channels()

    @property
    def named_channels(self):
        if self.__named_channels is None and self.channels is not None:
            self.__named_channels = {}
            for c in self.channels:
                self.__named_channels[c.name] = c

        return self.__named_channels

    @property
    def is_compressed(self):
        return self.graph_header.compressed

    @property
    def data_length(self):
        if self.is_compressed:
            return 0
        return sum([c.data_length for c in self.channels])

    def __str__(self):
        return("Biopac file (rev %s): %s channels, %s samples/sec" % (
            self.graph_header.file_revision, self.graph_header.channel_count,
            self.samples_per_second))

    def __repr__(self):
        return str(self)

    @property
    def time_index(self):
        if self.__time_index is not None:
            return self.__time_index

        total_samples = max(
            [ch.frequency_divider * ch.point_count
                for ch in self.channel_headers])
        total_seconds = total_samples / self.samples_per_second
        self.__time_index = np.linspace(0, total_seconds, total_samples)
        return self.__time_index

    def __build_channels(self):
        return [
            Channel.from_headers(
                ch, cdh, self.samples_per_second, self.time_index)
            for ch, cdh in zip(
                self.channel_headers, self.channel_dtype_headers)
        ]

    def __set_channel_time_indexes(self):
        for c in self.channels:
            c.time_index = self.time_index[::c.frequency_divider]


class Channel(object):
    """
    An individual channel of Biopac data. Has methods to access raw data from
    the file, as well as a scaled copy if the raw data is in integer format.
    Also generally created by reader.Reader.
    """

    def __init__(
            self,
            frequency_divider=None, raw_scale_factor=None, raw_offset=None,
            name=None, units=None, fmt_str=None, samples_per_second=None,
            point_count=None, time_index=None, order_num=None):

        self.frequency_divider = frequency_divider
        self.raw_scale_factor = raw_scale_factor
        self.raw_offset = raw_offset
        self.name = name
        self.units = units
        self.fmt_str = fmt_str
        self.samples_per_second = samples_per_second
        self.point_count = point_count
        self.dtype = np.dtype(fmt_str)
        self.time_index = time_index
        self.order_num = order_num

        # Don't allocate storage automatically -- this means we can read
        # only some channels or stream the data without putting all the data
        # in memory.
        self.raw_data = None
        self.__data = None
        self.__upsampled_data = None

    @classmethod
    def from_headers(
            cls,
            chan_hdr,
            dtype_hdr,
            samples_per_second,
            base_time_index=None):
        divider = chan_hdr.frequency_divider
        chan_samp_per_sec = samples_per_second / divider
        time_index = None
        if base_time_index is not None:
            time_index = base_time_index[::divider][0:chan_hdr.point_count]

        return cls(
            frequency_divider=divider,
            raw_scale_factor=chan_hdr.raw_scale,
            raw_offset=chan_hdr.raw_offset,
            name=chan_hdr.name,
            units=chan_hdr.units,
            fmt_str=dtype_hdr.numpy_dtype,
            samples_per_second=chan_samp_per_sec,
            point_count=chan_hdr.point_count,
            time_index=time_index,
            order_num=chan_hdr.order_num
        )

    def _allocate_raw_data(self):
        self.raw_data = np.zeros(self.point_count, dtype=self.dtype)

    @property
    def sample_size(self):
        """
        The size, in bytes, of one sample's worth of data.
        """
        return self.dtype.itemsize

    @property
    def data_length(self):
        """
        The size, in bytes, of the entire channel's raw data stream.
        """
        return self.sample_size * self.point_count

    @property
    def loaded(self):
        return self.raw_data is not None

    @property
    def data(self):
        """
        The channel's data, scaled by the raw_scale_factor and offset. These
        will be the values reported by AcqKnowledge. Note: only integer data
        types are scaled and offset.
        """
        if not self.loaded:
            return None
        if self.__data is not None:
            return self.__data
        scale_factor = self.raw_scale_factor
        raw_offset = self.raw_offset
        if self.dtype.kind == "i":
            self.__data = (self.raw_data * scale_factor) + raw_offset
            scale_factor = 1.0
            raw_offset = 0.0
        else:
            self.__data = self.raw_data
        return self.__data

    @property
    def upsampled_data(self):
        """
        The channel's data, sampled at the native frequency of the file.
        All channels will have the same number of points using this method,
        unless recording stopped in the middle of a block.
        Nearest-neighbor sampling is used.
        """
        if self.__upsampled_data is None:
            total_samples = self.data.shape[0]*self.frequency_divider
            self.__upsampled_data = self.data[
                np.arange(total_samples)//self.frequency_divider]
        return self.__upsampled_data

    def free_data(self):
        self.raw_data = None
        self.__data = None
        self.__upsampled_data = None

    def __str__(self):
        return("Channel %s: %s samples, %s samples/sec, loaded: %s" % (
            self.name,
            self.point_count,
            self.samples_per_second,
            self.loaded))

    def __repr__(self):
        return str(self)


class Marker(object):
    """
    A marker -- some kind of annotation for an AcqKnowledge file. They all
    have a sample index and some text, and more modern ones can be one of
    several styles (say, a flag or a star or a waveform start) and can be
    attached to a particular channel.
    """

    def __init__(
            self,
            sample_index,
            text,
            channel,
            style=None):

        self.sample_index = sample_index
        self.text = text
        self.channel = channel
        self.style = style
        super(Marker, self).__init__()

    def __eq__(self, other):
        return all([
            self.sample_index == other.sample_index,
            self.text == other.text,
            self.channel == other.channel,
            self.style == other.style
        ])

    def __str__(self):
        return("Marker {0}: sample index: {1} channel: {2} style: {3}".format(
            self.text, self.sample_index, self.channel, self.style))

    def __repr__(self):
        return str(self)
