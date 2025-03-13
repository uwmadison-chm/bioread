# coding: utf8
# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2025 Board of Regents of the University of Wisconsin System
#
# Written Nate Vack <njvack@wisc.edu> with research from John Ollinger
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison
#
# This is a very simple thing to write an AcqKnowledge file as tab-delimited
# text. It'll produce very large output, so use it with care.

import csv


def write_text(datafile, out_stream, channel_indexes, missing_val):
    writer = csv.writer(out_stream, delimiter=str("\t"))
    if not channel_indexes:
        channel_indexes = range(len(datafile.channels))
    chans = [datafile.channels[i] for i in channel_indexes]
    headers = ["time (s)"] + [f"{c.name} ({c.units})" for c in chans]
    writer.writerow(headers)

    for i, t in enumerate(datafile.time_index):
        rd = [t] + [data_or_blank(c, i, missing_val) for c in chans]
        writer.writerow(rd)


def data_or_blank(channel, index, missing_val):
    ci = index // channel.frequency_divider
    if index % channel.frequency_divider == 0 and ci < channel.point_count:
        return channel.data[ci]
    return missing_val
