# coding: utf8
# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2025 Board of Regents of the University of Wisconsin System
#
# Written Nate Vack <njvack@wisc.edu> with research from John Ollinger
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison

"""
Turns an .acq file into a .mat file with scipy.io.savemat
"""

import numpy as np


class MatlabWriter:
    def __init__(
        self, data=None, filename=None, compress=False, oned_as="row", data_only=False, single=False
    ):
        self.data = data
        self.filename = filename
        self.compress = compress
        self.oned_as = oned_as
        self.data_only = data_only
        self.single = single
        self.write_meta = not self.data_only

    @classmethod
    def write_file(cls, data, filename, compress=False, oned_as="row", data_only=False, single=False):
        writer = cls(data, filename, compress, oned_as, data_only, single)
        writer.write()

    def write(self):
        from scipy.io import savemat

        d = self.__build_dict(self.data)
        savemat(
            self.filename,
            d,
            format="5",
            do_compression=self.compress,
            oned_as=self.oned_as,
        )

    def __build_dict(self, data):
        d = {}
        d["samples_per_second"] = data.samples_per_second
        nc = len(data.channels)
        channels = np.zeros(nc, dtype="O")
        for i in range(nc):
            data_format = "=f4" if self.single else "=f8"
            c = data.channels[i]
            chan_dict = {}
            chan_dict["data"] = c.data.astype(data_format)
            chan_dict["samples_per_second"] = c.samples_per_second
            chan_dict["name"] = c.name
            chan_dict["frequency_divider"] = c.frequency_divider
            chan_dict["units"] = c.units
            channels[i] = chan_dict

        d["channels"] = channels

        if self.write_meta:
            d["event_markers"] = self.__build_markers(data)
            if data.journal is not None:
                d["journal"] = data.journal

        return d

    def __build_markers(self, data):
        markers = np.zeros(len(data.event_markers), dtype="O")
        for i, marker in enumerate(data.event_markers):
            md = {
                "label": marker.text,
                "sample_index": marker.sample_index,
                "type_code": marker.type_code or "None",
                "type": marker.type or "None",
                "date_created": marker.date_created_str,
                "channel_number": marker.channel_number or -1,
                "channel": marker.channel_name or "Global",
            }
            markers[i] = md
        return markers
