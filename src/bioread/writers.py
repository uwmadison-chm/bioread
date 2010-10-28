# coding: utf8
# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2010 Board of Regents of the University of Wisconsin System
#
# Written by John Ollinger <ollinger@wisc.edu> and Nate Vack <njvack@wisc.edu>
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison
# Project home: http://github.com/njvack/bioread

from scipy.io import savemat


class MatlabWriter(object):

    def __init__(self,
        data=None, filename=None, compress=False, oned_as='row'):
        self.data = data
        self.filename = filename
        self.compress = compress
        self.oned_as = oned_as

    @classmethod
    def write_file(cls, data, filename, compress=False,
        oned_as='row'):
        writer = cls(data, filename, compress, oned_as)
        writer.write()

    def write(self):
        d = self.__build_dict(self.data)
        savemat(self.filename, d,
            format='5', do_compression=self.compress,
            oned_as=self.oned_as)

    def __build_dict(self, data):
        d = {}
        d['samples_per_second'] = data.samples_per_second
        channels = {}
        channel_headers = {}
        channel_dtype_headers = {}
        for i in range(len(data.channels)):
            c = data.channels[i]
            label = "n%s" % (i+1)
            chan_dict = {}
            chan_dict['data'] = c.data
            chan_dict['samples_per_second'] = c.samples_per_second
            chan_dict['name'] = c.name
            chan_dict['frequency_divider'] = c.freq_divider
            chan_dict['units'] = c.units
            channels[label] = chan_dict
            channel_headers[label] = data.channel_headers[i].data
            channel_dtype_headers[label] = data.channel_dtype_headers[i].data

        d['channels'] = channels

        d['headers'] = {
            'graph': data.graph_header.data,
            'foreign': data.foreign_header.data,
            'channel': channel_headers,
            'channel_dtype': channel_dtype_headers}

        return d
