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
    
    
    
    # A list of tuples. Each tuple has three elements:
    # field_name, field_structure, min_version
    @property
    def __header_structure(self):
        return [
        ('nItemHeaderLen',          'h', 30),
        ('iVersion',                'i', 30),
        ]
    
    def __init__(self, acq_file):
        self.acq_file = acq_file
    
    def read(self):
        pass
    
    def _set_order_and_version(self):
        # Try unpacking the version string in both a bid and little-endian
        # fashion. Version string should be a small, positive integer.
        ver_hdr_part = self.__header_structure[0:2]
        ver_fmt_str = ''.join(e[1] for e in ver_hdr_part)
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
        
