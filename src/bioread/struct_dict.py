# coding: utf8
# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2010 Board of Regents of the University of Wisconsin System
#
# Written by John Ollinger <ollinger@wisc.edu> and Nate Vack <njvack@wisc.edu>
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison
# 

import struct


class StructDict(object):
    """
    This class allows you to declare a header's structure with name and size
    information, and then will unpack the struct into a dict. For example:

    >>> sd = StructDict('>', [('version', 'h'),('xy_dim', '2b'), ('name', '5s')])
    >>> sd.unpack('\x00\x01\x05\x10foo\x00\x00')
    {
        'version' : 1,
        'xy_dim' : (5, 16),
        'name' : 'foo'
    }
    """
    
    def __init__(self, byte_order_flag, struct_info=None):
        self.byte_order_flag = byte_order_flag
        self.struct_info = struct_info
        self.full_struct_info = None
    
    def unpack(self, data):
        self.__setup()
        unpacked = struct.unpack(self.format_string, data)
        output = {}
        for name, fs, start_index, end_index in self.full_struct_info:
            l = end_index-start_index
            if l == 1:
                val = unpacked[start_index]
            else:
                val = unpacked[start_index:end_index]
            output[name] = val
        return output
    
    @property
    def format_string(self):
        s = ''.join([si[1] for si in self.struct_info])
        return self.__bof_fs(s)

    def __setup(self):
        if self.full_struct_info is None:
            self.full_struct_info = self.__full_struct_info()
    
    def __bof_fs(self, format_str):
        return self.byte_order_flag + format_str
    
    def __unpacked_element_count(self, format_str):
        # We need to figure this out by actually faking some data and
        # unpacking it. Crazy, huh?
        f_str = self.__bof_fs(format_str)
        f_len = struct.calcsize(f_str)
        dummy = '\x00'*f_len
        unpacked = struct.unpack(f_str, dummy)
        return len(unpacked) # The number of elements in the tuple
    
    def __full_struct_info(self):
        full_struct_info = []
        start_index = 0
        end_index = 0
        for name, fs in self.struct_info:
            tup_len = self.__unpacked_element_count(fs)
            end_index = start_index + tup_len
            full_struct_info.append((name, fs, start_index, end_index))
            start_index = end_index
        return full_struct_info