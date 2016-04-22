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
        self.event_markers = None
        self.journal_header = None
        self.journal = None
        self.__named_channels = None
        self.__time_index = None

        self.channels = self.__build_channels()
        self.channel_order_map = dict(
            [[c.order_num, c] for c in self.channels]
        )

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
        return("AcqKnowledge file (rev %s): %s channels, %s samples/sec" % (
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
                chan_hdr=ch,
                dtype_hdr=cdh,
                samples_per_second=self.samples_per_second,
                datafile=self)
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
            point_count=None, order_num=None, datafile=None):

        self.frequency_divider = frequency_divider
        self.raw_scale_factor = raw_scale_factor
        self.raw_offset = raw_offset
        self.name = name
        self.units = units
        self.fmt_str = fmt_str
        self.samples_per_second = samples_per_second
        self.point_count = point_count
        self.dtype = np.dtype(fmt_str)
        # For some reason, scale and offset lie for float files
        if self.dtype.kind == 'f':
            self.raw_scale_factor = 1
            self.raw_offset = 0
        self.order_num = order_num
        self.datafile = datafile

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
            datafile=None):
        divider = chan_hdr.frequency_divider
        chan_samp_per_sec = samples_per_second / divider

        return cls(
            frequency_divider=divider,
            raw_scale_factor=chan_hdr.raw_scale,
            raw_offset=chan_hdr.raw_offset,
            name=chan_hdr.name,
            units=chan_hdr.units,
            fmt_str=dtype_hdr.numpy_dtype,
            samples_per_second=chan_samp_per_sec,
            point_count=chan_hdr.point_count,
            order_num=chan_hdr.order_num,
            datafile=datafile
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
    def time_index(self):
        if self.datafile is None:
            return None
        div = self.frequency_divider
        return self.datafile.time_index[::div][0:self.point_count]

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
        if self.dtype.kind == 'f':
            self.__data = self.raw_data
        else:
            self.__data = (
                (self.raw_data * self.raw_scale_factor) + self.raw_offset)
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


class EventMarker(object):
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
            channel_number,
            channel=None,
            type_code=None):

        self.sample_index = sample_index
        self.text = text
        self.channel_number = channel_number
        self.channel = channel
        self.type_code = type_code
        super(EventMarker, self).__init__()

    def __eq__(self, other):
        return all([
            self.sample_index == other.sample_index,
            self.text == other.text,
            self.channel_number == other.channel_number,
            self.type_code == other.type_code
        ])

    def __str__(self):
        return("EventMarker {}: idx: {} channel: {} type_code: {}".format(
            self.text,
            self.sample_index,
            self.channel_number,
            self.type_code))

    def __repr__(self):
        return str(self)

    @property
    def type(self):
        if self.type_code is None:
            return "None"
        return MARKER_TYPE_MAP.get(self.type_code, "Unknown")

    @property
    def channel_sample_index(self):
        if self.channel is None:
            return None
        return self.sample_index // self.channel.frequency_divider

    @property
    def channel_name(self):
        if self.channel is None:
            return None
        return self.channel.name


MARKER_TYPE_MAP = {
    'apnd': 'Append',
    'defl': 'Default',
    'wfon': 'Waveform Onset',
    'wfof': 'Waveform End',
    'nois': 'Change in Signal Quality',
    'rhyt': 'Change in Rhythm',
    'recv': 'Recovery',
    'max ': 'Maximum',
    'min ': 'Minimum',
    'rset': 'Reset',
    'cmlb': 'Communication Lost Begin',
    'cmle': 'Communication Lost End',
    'ansh': 'Short Arrow',
    'anmd': 'Medium Arrow',
    'anlg': 'Long Arrow',
    'flag': 'Flag',
    'star': 'Star',
    'usr1': 'User Type 1',
    'usr2': 'User Type 2',
    'usr3': 'User Type 3',
    'usr4': 'User Type 4',
    'usr5': 'User Type 5',
    'usr6': 'User Type 6',
    'usr7': 'User Type 7',
    'usr8': 'User Type 8',
    'usr9': 'User Type 9',
    'qrsb': 'QRS Onset',
    'qrs ': 'QRS Peak',
    'qrse': 'QRS End',
    'tbeg': 'T-wave Onset',
    't   ': 'T-wave Peak',
    'tend': 'T-wave End',
    'pbeg': 'P-wave Onset',
    'p   ': 'P-wave Peak',
    'pend': 'P-wave End',
    'q   ': 'Q-wave Peak',
    's   ': 'S-wave Peak',
    'u   ': 'U-wave Peak',
    'pq  ': 'PQ Junction',
    'jpt ': 'J-point',
    'stch': 'ST Segment Change',
    'tch ': 'T-wave Change',
    'nrml': 'Normal Beat',
    'pace': 'Paced Beat',
    'pfus': 'Fusion of Paced and Normal Beat',
    'lbbb': 'Left Bundle Branch Block Beat',
    'rbbb': 'Right Bundle Branch Block Beat',
    'bbb ': 'Bundle Branch Block Beat',
    'apc ': 'Atrial Premature Beat',
    'aber': 'Aberrated Atrial Prematuire Beat',
    'npc ': 'Nodal Premature Beat',
    'svpb': 'Supraventricular Premature Beat',
    'pvc ': 'Premature Ventricular Contraction',
    'ront': 'R-on-T Premature Ventricular Contraction',
    'fusi': 'Fusion of Ventricular and Normal Beat',
    'aesc': 'Atrial Escape Beat',
    'nesc': 'Nodal Escape Beat',
    'sves': 'Supraventricular Escape Beat',
    'vesc': 'Ventricular Escape Beat',
    'syst': 'Systole',
    'dias': 'Diastole',
    'edp ': 'End Diastolic Pressure',
    'aptz': 'A-point',
    'bptz': 'B-point',
    'cptz': 'C-point',
    'xptz': 'X-point',
    'yptz': 'Y-point',
    'optz': 'O-point',
    'plat': 'Plateau',
    'upst': 'Upstroke',
    'vfon': 'Start of Ventricular Flutter',
    'flwa': 'Ventricular Flutter Wave',
    'vfof': 'End of Ventricular Flutter',
    'pesp': 'Pacemaker Artifact',
    'arfc': 'Isolated QRS-like Artifact',
    'napc': 'Non-conducted P-wave',
    'base': 'Baseline',
    'dose': 'Dose',
    'wash': 'Wash',
    'apon': 'Spike Episode Begin',
    'apof': 'Spike Episode End',
    'rein': 'Inspire Start',
    'reot': 'Expire Start',
    'reap': 'Apnea Start',
    'stim': 'Stimulus Delivery',
    'resp': 'Response',
    'scr ': 'Skin Conductance Response',
    'sscr': 'Specific SCR',
    'ctr1': 'Cluster 1',
    'ctr2': 'Cluster 2',
    'ctr3': 'Cluster 3',
    'ctr4': 'Cluster 4',
    'ctr5': 'Cluster 5',
    'ctr6': 'Cluster 6',
    'ctr7': 'Cluster 7',
    'ctr8': 'Cluster 8',
    'ctr9': 'Cluster 9',
    'ctrn': 'Cluster n',
    'cend': 'End Cluster',
    'outl': 'Outlier',
    'tran': 'Training Set',
    'cut ': 'Cut',
    'vb  ': 'Paste Begin',
    've  ': 'Paste End',
    'selb': 'Selection Begin',
    'sele': 'Selection End',
    'steb': 'Start of Eye Blink Artifact',
    'eneb': 'End of Eye Blink Artifact',
    'sexc': 'Start of Excursion Artifact',
    'eexc': 'End of Excursion Artifact',
    'ssat': 'Start of Saturation Artifact',
    'esat': 'End of Saturation Artifact',
    'sspk': 'Start of Spike Artifact',
    'espk': 'End of Spike Artifact',
    'semg': 'Start of EMG Artifact',
    'eemg': 'End of EMG Artifact',
    'wles': 'Workload - EMG Start',
    'wlee': 'Workload - EMG End',
    'ipss': 'Workload - Invalid PSD Start',
    'ipse': 'Workload - Invalid PSD End',
    'ddst': 'Dummy Data Start',
    'dded': 'Dummy Data End',
    'idst': 'Misaligned Data',
    'bprs': 'Button Pressed',
    'leho': 'Left Eye Hit Object',
    'reho': 'Right Eye Hit Object',
    'smis': 'SMI Stimulus Image Has Been Presented to the Subject',
    'mors': 'Start Out of Range',
    'more': 'End Out of Range'
}
