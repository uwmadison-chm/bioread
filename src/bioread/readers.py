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
from file_versions import *

# TODO: Factor out a reader for Pre-Version-4 files. Or something?

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
    def _graph_header_structure(self):
        return [
        ('nItemHeaderLen'           ,'h'    ,V_ALL ),
        ('lVersion'                 ,'l'    ,V_ALL ),
        ('lExtItemHeaderLen'        ,'l'    ,V_20a ),
        ('nChannels'                ,'h'    ,V_20a ),
        ('nHorizAxisType'           ,'h'    ,V_20a ),
        ('nCurChannel'              ,'h'    ,V_20a ),
        ('dSampleTime'              ,'d'    ,V_20a ),
        ('dTimeOffset'              ,'d'    ,V_20a ),
        ('dTimeScale'               ,'d'    ,V_20a ),
        ('dTimeCursor1'             ,'d'    ,V_20a ),
        ('dTimeCursor2'             ,'d'    ,V_20a ),
        ('rcWindow'                 ,'4h'   ,V_20a ),
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
        ('BShowToolBar'             ,'h'    ,V_30r ),
        ('BShowChannelButtons'      ,'h'    ,V_30r ),
        ('BShowMeasurements'        ,'h'    ,V_30r ),
        ('BShowMarkers'             ,'h'    ,V_30r ),
        ('BShowJournal'             ,'h'    ,V_30r ),
        ('CurXChannel'              ,'h'    ,V_30r ),
        ('MmtPrecision'             ,'h'    ,V_30r ),
        ('NMeasurementRows'         ,'h'    ,V_303 ),
        ('mmt40'                    ,'40h'  ,V_303 ),
        ('mmtChan40'                ,'40h'  ,V_303 ),
        ('MmtCalcOpnd1'             ,'40h'  ,V_35x ),
        ('MmtCalcOpnd2'             ,'40h'  ,V_35x ),
        ('MmtCalcOp'                ,'40h'  ,V_35x ),
        ('MmtCalcConstant'          ,'40d'  ,V_35x ),
        ('bNewGridWithMinor'        ,'l'    ,V_370 ),
        ('colorMajorGrid'           ,'4B'   ,V_370 ),
        ('colorMinorGrid'           ,'4B'   ,V_370 ),
        ('wMajorGridStyle'          ,'h'    ,V_370 ),
        ('wMinorGridStyle'          ,'h'    ,V_370 ),
        ('wMajorGridWidth'          ,'h'    ,V_370 ),
        ('wMinorGridWidth'          ,'h'    ,V_370 ),
        ('bFixedUnitsDiv'           ,'l'    ,V_370 ),
        ('bMid_Range_Show'          ,'l'    ,V_370 ),
        ('dStart_Middle_Point'      ,'d'    ,V_370 ),
        ('dOffset_Point'            ,'d'    ,V_370 ),
        ('hGrid'                    ,'d'    ,V_370 ),
        ('vGrid'                    ,'d'    ,V_370 ),
        ('bEnableWaveTools'         ,'l'    ,V_370 ),
        ('Reserved'                 ,'20b'  ,V_381 ),
        ('bOverlapMode'             ,'l'    ,V_381 ),
        ('bShowHardware'            ,'l'    ,V_381 ),
        ('bXAutoPlot'               ,'l'    ,V_381 ),
        ('bXAutoScroll'             ,'l'    ,V_381 ),
        ('bStartButtonVisible'      ,'l'    ,V_381 ),
        ('bCompressed'              ,'l'    ,V_381 ),
        ('bAlwaysStartButtonVisible','l'    ,V_381 ),
        ('pathVideo'                ,'260s' ,V_382 ),
        ('optSyncDelay'             ,'l'    ,V_382 ),
        ('syncDelay'                ,'d'    ,V_382 ),
        ('bHRP_PasteMeasurements'   ,'l'    ,V_382 ),
        ('graphType'                ,'l'    ,V_390 ),
        ('mmtCalcExpr'              ,'10240s',V_390 ),
        ('mmtMomentOrder'           ,'40l'  ,V_390 ),
        ('mmtTimeDelay'             ,'40l'  ,V_390 ),
        ('mmtEmbedDim'              ,'40l'  ,V_390 ),
        ('mmtMIDelay'               ,'40l'  ,V_390 ),
        ]
    
    @property
    def _channel_header_structure(self):
        return [
        ('lChanHeaderLen'           ,'l'    ,V_20a ),
        ('nNum'                     ,'h'    ,V_20a ),
        ('szCommentText'            ,'40s'  ,V_20a ),
        ('rgbColor'                 ,'4B'   ,V_20a ),
        ('nDispChan'                ,'h'    ,V_20a ),
        ('dVoltOffset'              ,'d'    ,V_20a ),
        ('dVoltScale'               ,'d'    ,V_20a ),
        ('szUnitsText'              ,'20s'  ,V_20a ),
        ('lBufLength'               ,'l'    ,V_20a ),
        ('dAmplScale'               ,'d'    ,V_20a ),
        ('dAmplOffset'              ,'d'    ,V_20a ),
        ('nChanOrder'               ,'h'    ,V_20a ),
        ('nDispSize'                ,'h'    ,V_20a )
        ]
    
    def _headers_for(self, hstruct, ver):
        return [hs for hs in hstruct if hs[2] <= ver]
            
    def _graph_headers_for(self, ver):
        return self._headers_for(self._graph_header_structure, ver)
    
    def read(self):
        self.__setup()
        self.graph_header = self.__read_header(self.graph_header_sd, 0)
        self.channel_headers = []
        channel_head_len = 0 # This will get changed!
        for i in range(self.graph_header['nChannels']):
            channel_offset = (
                self.graph_header['lExtItemHeaderLen'] + i*channel_head_len)
            channel_header = self.__read_header(
                self.channel_header_sd, channel_offset)
            # This will be the same for all channels; easier to set it for all
            channel_head_len = channel_header['lChanHeaderLen']
            self.channel_headers.append(channel_header)
        
    def __setup(self):
        if self.byte_order_flag is not None:
            return
        # TODO: Extract this into a factory class
        self.__set_order_and_version()
        self.__build_header_struct_dicts()
    
    def __build_header_struct_dicts(self):
        graph_header = self._headers_for(
            self._graph_header_structure, self.file_version)
        self.graph_header_sd = StructDict(self.byte_order_flag, graph_header)

        channel_header = self._headers_for(
            self._channel_header_structure, self.file_version)
        self.channel_header_sd = StructDict(
            self.byte_order_flag, channel_header)
    
    def __read_header(self, struct_dict, offset):
        self.acq_file.seek(offset)
        data = self.acq_file.read(struct_dict.len_bytes)
        return struct_dict.unpack(data)
        
    def __set_order_and_version(self):
        # Try unpacking the version string in both a bid and little-endian
        # fashion. Version string should be a small, positive integer.
        self.acq_file.seek(0)
        mh = self._graph_headers_for(V_ALL)
        ver_fmt_str = ''.join([s[1] for s in mh])
        ver_len = struct.calcsize('<'+ver_fmt_str)
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