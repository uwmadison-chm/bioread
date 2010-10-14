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
    # field_name, field_structure, first_version_included
    @property
    def _main_header_structure(self):
        from file_versions import V_ALL, V_20a
        return [
        ('nItemHeaderLen'           ,'h'    ,V_ALL ),
        ('iVersion'                 ,'i'    ,V_ALL ),
        ('iExtItemHeaderLen'        ,'i'    ,V_20a ),
        ('nChannels'                ,'h'    ,V_20a ),
        ('nHorizAxisType'           ,'h'    ,V_20a ),
        ('nCurChannel'              ,'h'    ,V_20a ),
        ('dSampleTime'              ,'d'    ,V_20a ),
        ('dTimeOffset'              ,'d'    ,V_20a ),
        ('dTimeScale'               ,'d'    ,V_20a ),
        ('dTimeCursor1'             ,'d'    ,V_20a ),
        ('dTimeCursor2'             ,'d'    ,V_20a ),
        ('rcWindow'                 ,'8b'   ,V_20a ),
        ('nMeasurement'             ,'6h'   ,V_20a ),
        ('fHilite'                  ,'h'    ,V_20a ),
        ('dFirstTimeOffset'         ,'d'    ,V_20a ),
        ('nRescale'                 ,'h'    ,V_20a ),
        ('szHorizUnits1'            ,'40s'  ,V_20a ),
        ('szHorizUnits2'            ,'10s'  ,V_20a ),
        ('nInMemory'                ,'h'    ,V_20a ),
        ('fGrid'                    ,'h'    ,V_20a ),
        ('fMarkers'                 ,'h'    ,V_20a ),
        ('nPlotDraft'               ,'h'    ,V_20a ),
        ('nDispMode'                ,'h'    ,V_20a ),
        ('rRReserved'               ,'h'    ,V_20a ),
        ]
    
    def _headers_for(self, hstruct, ver):
        return [hs for hs in hstruct if hs[2] >= ver]
            
    def _main_headers_for(self, ver):
        return self._headers_for(self._main_header_structure, ver)
    
    def read(self):
        pass
    
    def _set_order_and_version(self):
        # Try unpacking the version string in both a bid and little-endian
        # fashion. Version string should be a small, positive integer.
        ver_fmt_str = self._main_header_str_for(0, False)
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
        

def header_name_offsets(header_info):
    """ 
    For debugging purposes -- the documenataion from Biopac reports
    offsets. This will let us know when we're following spec.
    """
    tab = []
    for i in range(len(header_info)):
        hparts = header_info[0:i]
        fmt_str = '<'+''.join([hp[1] for hp in hparts])
        tab.append((header_info[i][0], struct.calcsize(fmt_str)))
    return tab
        