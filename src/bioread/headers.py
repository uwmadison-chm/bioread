# coding: utf8
# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2010 Board of Regents of the University of Wisconsin System
#
# Written by John Ollinger <ollinger@wisc.edu> and Nate Vack <njvack@wisc.edu>
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison

import struct

from struct_dict import StructDict


class Header(object):
    """
    Represents a binary header and its unpacking
    """
    
    def __init__(self, struct_dict):
        self.struct_dict = struct_dict
        self.offset = None
        self.unpacked = {}
    
    def unpack_from_str(self, str_data):
        self.raw_data = str_data
        return self.__unpack_data()
        
    def unpack_from_file(self, data_file, offset):
        self.offset = offset
        data_file.seek(offset)
        self.raw_data = data_file.read(self.struct_dict.len_bytes)
        return self.__unpack_data()
    
    def __unpack_data(self):
        self.data = self.struct_dict.unpack(self.raw_data)
        return self.data
    
    def __getitem__(self, key):
        return self.data[key]
    
    def __str__(self):
        return str(self.data)

