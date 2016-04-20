# coding: utf8
# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2016 Board of Regents of the University of Wisconsin System
#
# Written Nate Vack <njvack@wisc.edu> with research from John Ollinger
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison
# Project home: http://github.com/njvack/bioread
#
# This is a very simple thing to write an AcqKnowledge file as tab-delimited
# text. It'll produce very large output, so use it with care.

from __future__ import unicode_literals
import csv


def write_text(datafile, out_stream, channel_indexes, missing_val):
    writer = csv.writer(out_stream, delimiter=str("\t"))
    if not channel_indexes:
        channel_indexes = range(len(datafile.channels))
    chans = [datafile.channels[i] for i in channel_indexes]
    headers = ["time (s)"] + [
        "{0} ({1})".format(c.name, c.units) for c in chans]
    headers = [s.encode('utf-8') for s in headers]
    writer.writerow(headers)
    for i, t in enumerate(datafile.time_index):
        rd = [t] + [data_or_blank(c, i, missing_val) for c in chans]
        writer.writerow(rd)


def data_or_blank(channel, index, missing_val):
    ci = index // channel.frequency_divider
    if index % channel.frequency_divider == 0 and ci < channel.point_count:
        return channel.data[ci]
    return missing_val
