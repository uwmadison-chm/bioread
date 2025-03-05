# coding: utf8
# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2025 Board of Regents of the University of Wisconsin System
#
# Written Nate Vack <njvack@wisc.edu> with research from John Ollinger
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison
# Project home: http://github.com/njvack/bioread

import numpy as np


class MatlabWriter(object):

    def __init__(
            self,
            data=None,
            filename=None,
            compress=False,
            oned_as='row',
            data_only=False):
        self.data = data
        self.filename = filename
        self.compress = compress
        self.oned_as = oned_as
        self.data_only = data_only

    @classmethod
    def write_file(
            cls, data, filename, compress=False, oned_as='row', data_only=False):
        writer = cls(data, filename, compress, oned_as, data_only)
        writer.write()

    def write(self):
        from scipy.io import savemat
        d = self.__build_dict(self.data)
        savemat(
            self.filename, d,
            format='5', do_compression=self.compress,
            oned_as=self.oned_as)

    def __build_dict(self, data):
        d = {}
        d['samples_per_second'] = data.samples_per_second
        nc = len(data.channels)
        channels = np.zeros(nc, dtype='O')
        channel_headers = np.zeros(nc, dtype='O')
        channel_dtype_headers = np.zeros(nc, dtype='O')
        for i in range(nc):
            c = data.channels[i]
            chan_dict = {}
            chan_dict['data'] = c.data.astype("=f8")
            chan_dict['samples_per_second'] = c.samples_per_second
            chan_dict['name'] = c.name
            chan_dict['frequency_divider'] = c.frequency_divider
            chan_dict['units'] = c.units
            channels[i] = chan_dict
            channel_headers[i] = data.channel_headers[i].data
            channel_dtype_headers[i] = data.channel_dtype_headers[i].data

        d['channels'] = channels

        if not self.data_only:
            d['event_markers'] = self.__build_markers(data)

        d['headers'] = {
            'graph': data.graph_header.data,
            'foreign': data.foreign_header.data,
            'channel': channel_headers,
            'channel_dtype': channel_dtype_headers}

        return d

    def __build_markers(self, data):
        markers = np.zeros(len(data.event_markers), dtype='O')
        for i, marker in enumerate(data.event_markers):
            md = {
                'label': marker.text,
                'sample_index': marker.sample_index,
                'type_code': marker.type_code or 'None',
                'type': marker.type or 'None',
                'date_created': marker.date_created_str,
                'channel_number': marker.channel_number or -1,
                'channel': marker.channel_name or 'Global'
            }
            markers[i] = md
        return markers
