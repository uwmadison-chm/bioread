# coding: utf8
# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2016 Board of Regents of the University of Wisconsin System
#
# Written Nate Vack <njvack@wisc.edu> with research from John Ollinger
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison
# Project home: http://github.com/njvack/bioread

from __future__ import absolute_import

from bioread.vendor.ordereddict import OrderedDict

import struct


class StructDict(object):
    """
    TODO: This class is basically a shittier version of ctypes. Remove this
          and replace it with ctypes.

    This class allows you to declare a binary file's header structure with
    name and type information, and then will unpack the header into a
    dictionary.
    For example:
    >>> header_structure = [
        ('version', 'h'), ('xy_dim', '2b'), ('name', '5s')
    ]
    >>> sd = StructDict('>', header_structure)
    >>> sd.format_string
    'h2b5s'
    >>> sd.len_bytes
    9
    >>> header_data = '\x00\x01\x05\x10foo\x00\x00'
    >>> sd.unpack(header_data)
    {
        'version' : 1,     # Single-length elements are de-tupelized
        'xy_dim' : (5, 16),
        'name' : 'foo'     # Strings get trailing nulls trimmed
    }
    """

    def __init__(self, byte_order_char, struct_info=None):
        self.byte_order_char = byte_order_char
        self.struct_info = struct_info
        self.full_struct_info = None

    def unpack(self, data):
        """
        Return a dict with the unpacked data.
        """
        self.__setup()
        unpacked = struct.unpack(self.format_string, data)
        output = OrderedDict()
        for name, fs, start_index, end_index in self.full_struct_info:
            l = end_index-start_index
            if l == 1:
                val = unpacked[start_index]
            else:
                val = unpacked[start_index:end_index]
            if type(val) == str:
                null_idx = val.find("\x00")
                if null_idx == -1:
                    null_idx = len(val)
                val = val[0:null_idx]
            output[name] = val
        return output

    def labeled_offsets_lengths(self):
        """
        Primarily for debugging purposes: generate a list of byte offsets
        and struct lengths, so you can see what fields are where. Then you
        can explore with your hex editor, or compare against spec, or
        whatever.

        Example:

        >>> sd = StructDict('>', [('version', 'h'), ('header_len', 'l')])
        >>> sd.labeled_offsete_lengths()
        [
            ('version', 'h', 0, 2),
            ('header_len', 'l', 2, 4)
        ]

        """
        table = []
        build_fs = self.byte_order_char
        for si in self.struct_info:
            name, fs = si[0:2]
            f_offset = struct.calcsize(build_fs)
            f_len = struct.calcsize(self.byte_order_char+fs)
            build_fs += fs
            table.append((name, fs, f_offset, f_len))

        return table

    @property
    def len_bytes(self):
        return struct.calcsize(self.format_string)

    @property
    def len_elements(self):
        return len(self.struct_info)

    @property
    def format_string(self):
        s = ''.join([si[1] for si in self.struct_info])
        return self.__bof_fs(s)

    def __setup(self):
        if self.full_struct_info is None:
            self.full_struct_info = self.__full_struct_info()

    def __bof_fs(self, format_str):
        return self.byte_order_char + format_str

    def __unpacked_element_count(self, format_str):
        # We need to figure this out by actually faking some data and
        # unpacking it. Crazy, huh?
        f_str = self.__bof_fs(format_str)
        f_len = struct.calcsize(f_str)
        dummy = struct.pack('%ss' % f_len, b'')
        unpacked = struct.unpack(f_str, dummy)
        return len(unpacked) # The number of elements in the tuple

    def __full_struct_info(self):
        full_struct_info = []
        start_index = 0
        end_index = 0
        for si in self.struct_info:
            name, fs = si[0:2]
            tup_len = self.__unpacked_element_count(fs)
            end_index = start_index + tup_len
            full_struct_info.append((name, fs, start_index, end_index))
            start_index = end_index
        return full_struct_info
