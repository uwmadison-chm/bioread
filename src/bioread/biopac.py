# coding: utf8
# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2010 Board of Regents of the University of Wisconsin System
#
# Written by John Ollinger <ollinger@wisc.edu> and Nate Vack <njvack@wisc.edu>
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison
# Project home: http://github.com/njvack/bioread

import numpy as np


class Datafile(object):
    """
    A data file for the AcqKnowledge system. Generally, gonna be created
    from a file by readers.AcqReader.
    """

    def __init__(self,
            graph_header=None, channel_headers=None, foreign_header=None,
            channel_dtype_headers=None, samples_per_second=None, name=None):
        self.graph_header = graph_header
        self.channel_headers = channel_headers
        self.foreign_header = foreign_header
        self.channel_dtype_headers = channel_dtype_headers
        self.samples_per_second = samples_per_second
        self.name = name
        self.channels = None
        self.__named_channels = None

    @property
    def named_channels(self):
        if self.__named_channels is None and self.channels is not None:
            self.__named_channels = {}
            for c in self.channels:
                self.__named_channels[c.name] = c

        return self.__named_channels

    def __unicode__(self):
        return("Biopac file (rev %s): %s channels, %s samples/sec" % (
            self.graph_header.file_revision, len(self.channels),
            self.samples_per_second))

    def __str__(self):
        return str(unicode(self))

    def __repr__(self):
        return str(self)


class Channel(object):
    """
    An individual channel of Biopac data. Has methods to access raw data from
    the file, as well as a scaled copy if the raw data is in integer format.
    Also generally created by readers.AcqReader.
    """

    def __init__(self,
        freq_divider=None, raw_scale_factor=None, raw_offset=None,
        raw_data=None, name=None, units=None, fmt_str=None,
        samples_per_second=None):

        self.freq_divider = freq_divider
        self.raw_scale_factor = raw_scale_factor
        self.raw_offset = raw_offset
        self.name = name
        self.units = units
        self.fmt_str = fmt_str
        self.samples_per_second = samples_per_second
        self.raw_data = raw_data
        self.__data = None
        self.__upsampled_data = None

    def should_sample_at(self, base_index):
        """
        Return true if the channel should be sampled in the graph file's
        base_index-th sample. This is when we're base_index is exactly
        divisible by freq_counter and we have room for it in our data.
        """
        return (
            (base_index % self.freq_divider) == 0 and
            (base_index // self.freq_divider) < self.point_count)

    @property
    def sample_size(self):
        """
        The size, in bytes, of one sample's worth of data.
        """
        return self.raw_data.dtype.itemsize

    @property
    def data_length(self):
        """
        The size, in bytes, of the entire channel's raw data stream.
        """
        return self.sample_size * self.point_count

    @property
    def data_proportion(self):
        """
        The sample size divided by the frequency divider.
        """
        return float(self.sample_size)/self.freq_divider

    @property
    def point_count(self):
        """
        Shorthand for len(self.raw_data).
        """
        return self.raw_data.shape[0]

    @property
    def data(self):
        """
        The channel's data, scaled by the raw_scale_factor and offset. These
        will be the values reported by AcqKnowledge. Note: only integer data
        types are scaled and offset.
        """
        scale_factor = self.raw_scale_factor
        raw_offset = self.raw_offset
        if self.fmt_str.find("d") >= 0: # test for float-ness
            scale_factor = 1.0
            raw_offset = 0.0
        if self.__data is None:
            self.__data = (self.raw_data*scale_factor)+raw_offset
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
            total_samples = self.data.shape[0]*self.freq_divider
            self.__upsampled_data = self.data[
                np.arange(total_samples)//self.freq_divider]
        return self.__upsampled_data

    def __unicode__(self):
        return("Channel %s: %s samples, %s samples/sec" % (
            self.name, len(self.raw_data), self.samples_per_second))

    def __str__(self):
        return str(unicode(self))

    def __repr__(self):
        return str(self)
