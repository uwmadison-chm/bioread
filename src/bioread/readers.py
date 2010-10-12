# coding: utf8
# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2010 Board of Regents of the University of Wisconsin System
#
# Written by John Ollinger <ollinger@wisc.edu> and Nate Vack <njvack@wisc.edu>
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison

import struct


class AcqReader(object):
    """ 
    Main class for reading AcqKnowledge files. You'll probably call it like:
    >>> data = AcqReader.read("some_file.acq")
    """
    
    def __init__(self, acq_file):
        self.acq_file = acq_file
        # This must be set by _set_order_and_version
        self.byte_order_flag = None
    
    # A list of tuples. Each tuple has three elements:
    # field_name, field_structure, min_version
    @property
    def __main_header_structure(self):
        return [
        ('nItemHeaderLen',          'h', 0),
        ('iVersion',                'i', 0),
        ]
    
    def __headers_for(self, hstruct, ver):
        return [hs for hs in hstruct if hs[2] >= ver]
        
    def __header_str_for(self, hstruct, ver, include_bof=True):
        """ Generates a format string for header data such as '>hi' """
        fmt_str = ''
        if include_bof:
            fmt_str = self.byte_order_flag
        fmt_str += ''.join([hs[1] for hs in self.__headers_for(hstruct, ver)])
        return fmt_str
    
    def __main_header_str_for(self, ver, include_bof=True):
        return self.__header_str_for(
            self.__main_header_structure, ver, include_bof)
    
    def __main_headers_for(self, ver):
        return self.__headers_for(self.__main_header_structure, ver)
    
    def read(self):
        pass
    
    def _set_order_and_version(self):
        # Try unpacking the version string in both a bid and little-endian
        # fashion. Version string should be a small, positive integer.
        ver_fmt_str = self.__main_header_str_for(0, False)
        ver_len = struct.calcsize('<'+ver_fmt_str)
        self.acq_file.seek(0)
        ver_data = self.acq_file.read(ver_len)
        
        byte_order_flags = ['<', '>']
        # Try both ways.
        byte_order_versions = [
            (struct.unpack(bof+ver_fmt_str, ver_data)[1], bof) for 
                bof in byte_order_flags
        ]
        
        # Limit to positive numbers, choose smallest.
        byte_order_versions = sorted([
            bp for bp in byte_order_versions if bp[0] > 0
        ])
        bp = byte_order_versions[0]
        
        self.byte_order_flag = bp[1]
        self.file_version = bp[0]
        
