# coding: utf8
# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2010 Board of Regents of the University of Wisconsin System
#
# Written by John Ollinger <ollinger@wisc.edu> and Nate Vack <njvack@wisc.edu>
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison

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
    
class Channel(object):
    """
    An individual channel of Biopac data. Has methods to access raw data from
    the file, as well as a scaled copy if the raw data is in integer format.
    Also generally created by readers.AcqReader.
    """
    def __init__(self, 
        freq_divider=None, raw_scale_factor=None, raw_offset=None, 
        raw_data=None, name=None, units=None, fmt_str=None):
        
        self.freq_divider = freq_divider
        self.raw_scale_factor = raw_scale_factor
        self.raw_offset = raw_offset
        self.raw_data = raw_data
        self.name = name
        self.units = units
        self.fmt_str = fmt_str
