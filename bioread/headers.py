# coding: utf8
# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2025 Board of Regents of the University of Wisconsin System
#
# Written Nate Vack <njvack@wisc.edu> with research from John Ollinger
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison
# Project home: http://github.com/njvack/bioread

# flake8: noqa: E203, E231, E741, E405, E403

# If I were to do this over again, I'd probably design it differently, without
# the somewhat weird versioning structure, as well as using ctypes instead
# of the weird StructDict thing we have.
# But, you know what? It's not worth changing this now, unless there turns out
# to be a compelling reason in the future.
# But if you're ever looking to do something similar, probably don't use this as
# a template of good code.


from bioread.struct_dict import StructDict
from bioread.file_revisions import *


class Header(object):
    """
    Represents a binary header and its unpacking
    """

    def __init__(self, struct_dict, encoding='utf-8'):
        self.struct_dict = struct_dict
        self.encoding = encoding
        self.offset = None
        self.unpacked = {}

    def unpack_from_str(self, str_data):
        self.raw_data = str_data
        self.__unpack_data()

    def unpack_from_file(self, data_file, offset):
        self.offset = offset
        data_file.seek(offset)
        self.raw_data = data_file.read(self.struct_dict.len_bytes)
        self.__unpack_data()

    @property
    def effective_len_bytes(self):
        """
        This will be overridden frequently -- it's used in navigating files.
        """
        return self.struct_dict.len_bytes

    @property
    def format_string(self):
        return self.struct_dict.format_string

    def __unpack_data(self):
        self.data = self.struct_dict.unpack(self.raw_data)

    def __getitem__(self, key):
        return self.data[key]

    def __str__(self):
        return str(self.data)


class VersionedHeaderStructure(object):

    def __init__(self, *structure_elements):
        self.structure_elements = structure_elements

    def elements_for(self, version):
        return [se for se in self.structure_elements if se[2] <= version]


class BiopacHeader(Header):
    """
    A simple superclass for GraphHeader, ChannelHeader, and friends.
    """

    def __init__(self, header_structure, file_revision, byte_order_char,
                 **kwargs):
        self.file_revision = file_revision
        self.byte_order_char = byte_order_char
        self.header_structure = header_structure
        sd = StructDict(byte_order_char,
                        header_structure.elements_for(file_revision))
        super().__init__(sd, **kwargs)


class GraphHeader(BiopacHeader):
    """
    The main Graph Header for an AcqKnowledge file. The documentation is much
    more complete for older files (revision 45 and below); I've looked to find
    the reliable (and essential) structure information in later files.

    The fields we really need to decode the data:
    Length
    Number of channels
    Sample time
    Data compressed?
    """

    def __init__(self, file_revision, byte_order_char, **kwargs):
        self.file_revision = file_revision
        super().__init__(self.__h_elts, file_revision, byte_order_char,
                         **kwargs)

    @property
    def effective_len_bytes(self):
        l = self.data['lExtItemHeaderLen']
        # We're going to skip an extra 40 bytes. This is kind of an ugly hack;
        # there seems to be a whole other 40-byte header in there.
        # if self.file_revision >= V_430:
        #     l += 40
        return l

    @property
    def channel_count(self):
        return self.data['nChannels']

    @property
    def sample_time(self):
        return self.data['dSampleTime']

    @property
    def compressed(self):
        return self.data.get('bCompressed') == 1

    @property
    def data_format(self):
        fmt = 'uncompressed'
        if self.compressed:
            fmt = 'compressed'
        return fmt

    @property
    def __version_bin(self):
        bin = 'Unknown'
        if self.file_revision < V_400B:
            bin = 'PRE_4'
        else:
            bin = 'POST_4'
        return bin

    @property
    def expected_padding_headers(self):
        if self.file_revision >= V_430:
            return self.data['hExpectedPaddings']
        return 0

    @property
    def __h_elt_versions(self):
        return {
            'PRE_4' : VersionedHeaderStructure(
                ('nItemHeaderLen'           ,'h'    ,V_ALL),
                ('lVersion'                 ,'l'    ,V_ALL),
                ('lExtItemHeaderLen'        ,'l'    ,V_20a),
                ('nChannels'                ,'h'    ,V_20a),
                ('nHorizAxisType'           ,'h'    ,V_20a),
                ('nCurChannel'              ,'h'    ,V_20a),
                ('dSampleTime'              ,'d'    ,V_20a),
                ('dTimeOffset'              ,'d'    ,V_20a),
                ('dTimeScale'               ,'d'    ,V_20a),
                ('dTimeCursor1'             ,'d'    ,V_20a),
                ('dTimeCursor2'             ,'d'    ,V_20a),
                ('rcWindow'                 ,'4h'   ,V_20a),
                ('nMeasurement'             ,'6h'   ,V_20a),
                ('fHilite'                  ,'h'    ,V_20a),
                ('dFirstTimeOffset'         ,'d'    ,V_20a),
                ('nRescale'                 ,'h'    ,V_20a),
                ('szHorizUnits1'            ,'40s'  ,V_20a),
                ('szHorizUnits2'            ,'10s'  ,V_20a),
                ('nInMemory'                ,'h'    ,V_20a),
                ('fGrid'                    ,'h'    ,V_20a),
                ('fMarkers'                 ,'h'    ,V_20a),
                ('nPlotDraft'               ,'h'    ,V_20a),
                ('nDispMode'                ,'h'    ,V_20a),
                ('rRReserved'               ,'h'    ,V_20a),
                ('BShowToolBar'             ,'h'    ,V_30r),
                ('BShowChannelButtons'      ,'h'    ,V_30r),
                ('BShowMeasurements'        ,'h'    ,V_30r),
                ('BShowMarkers'             ,'h'    ,V_30r),
                ('BShowJournal'             ,'h'    ,V_30r),
                ('CurXChannel'              ,'h'    ,V_30r),
                ('MmtPrecision'             ,'h'    ,V_30r),
                ('NMeasurementRows'         ,'h'    ,V_303),
                ('mmt40'                    ,'40h'  ,V_303),
                ('mmtChan40'                ,'40h'  ,V_303),
                ('MmtCalcOpnd1'             ,'40h'  ,V_35x),
                ('MmtCalcOpnd2'             ,'40h'  ,V_35x),
                ('MmtCalcOp'                ,'40h'  ,V_35x),
                ('MmtCalcConstant'          ,'40d'  ,V_35x),
                ('bNewGridWithMinor'        ,'l'    ,V_370),
                ('colorMajorGrid'           ,'4B'   ,V_370),
                ('colorMinorGrid'           ,'4B'   ,V_370),
                ('wMajorGridStyle'          ,'h'    ,V_370),
                ('wMinorGridStyle'          ,'h'    ,V_370),
                ('wMajorGridWidth'          ,'h'    ,V_370),
                ('wMinorGridWidth'          ,'h'    ,V_370),
                ('bFixedUnitsDiv'           ,'l'    ,V_370),
                ('bMid_Range_Show'          ,'l'    ,V_370),
                ('dStart_Middle_Point'      ,'d'    ,V_370),
                ('dOffset_Point'            ,'60d'  ,V_370),
                ('hGrid'                    ,'d'    ,V_370),
                ('vGrid'                    ,'60d'  ,V_370),
                ('bEnableWaveTools'         ,'l'    ,V_370),
                ('hozizPrecision'           ,'h'    ,V_373),
                ('Reserved'                 ,'20b'  ,V_381),
                ('bOverlapMode'             ,'l'    ,V_381),
                ('bShowHardware'            ,'l'    ,V_381),
                ('bXAutoPlot'               ,'l'    ,V_381),
                ('bXAutoScroll'             ,'l'    ,V_381),
                ('bStartButtonVisible'      ,'l'    ,V_381),
                ('bCompressed'              ,'l'    ,V_381),
                ('bAlwaysStartButtonVisible','l'    ,V_381),
                ('pathVideo'                ,'260s' ,V_382),
                ('optSyncDelay'             ,'l'    ,V_382),
                ('syncDelay'                ,'d'    ,V_382),
                ('bHRP_PasteMeasurements'   ,'l'    ,V_382),
                ('graphType'                ,'l'    ,V_390),
                ('mmtCalcExpr'              ,'10240s',V_390),
                ('mmtMomentOrder'           ,'40l'  ,V_390),
                ('mmtTimeDelay'             ,'40l'  ,V_390),
                ('mmtEmbedDim'              ,'40l'  ,V_390),
                ('mmtMIDelay'               ,'40l'  ,V_390),
            ),
            'POST_4' : VersionedHeaderStructure(
                ('nItemHeaderLen'           ,'h'    ,V_ALL),
                ('lVersion'                 ,'l'    ,V_ALL),
                ('lExtItemHeaderLen'        ,'l'    ,V_20a),
                ('nChannels'                ,'h'    ,V_20a),
                ('nHorizAxisType'           ,'h'    ,V_20a),
                ('nCurChannel'              ,'h'    ,V_20a),
                ('dSampleTime'              ,'d'    ,V_20a),
                ('dTimeOffset'              ,'d'    ,V_20a),
                ('dTimeScale'               ,'d'    ,V_20a),
                ('dTimeCursor1'             ,'d'    ,V_20a),
                ('dTimeCursor2'             ,'d'    ,V_20a),
                ('rcWindow'                 ,'4h'   ,V_20a),
                ('nMeasurement'             ,'6h'   ,V_20a),
                ('fHilite'                  ,'h'    ,V_20a),
                ('dFirstTimeOffset'         ,'d'    ,V_20a),
                ('nRescale'                 ,'h'    ,V_20a),
                ('szHorizUnits1'            ,'40s'  ,V_20a),
                ('szHorizUnits2'            ,'10s'  ,V_20a),
                ('nInMemory'                ,'h'    ,V_20a),
                ('fGrid'                    ,'h'    ,V_20a),
                ('fMarkers'                 ,'h'    ,V_20a),
                ('nPlotDraft'               ,'h'    ,V_20a),
                ('nDispMode'                ,'h'    ,V_20a),
                ('rRReserved'               ,'h'    ,V_20a),
                ('Unknown'                  ,'822B' ,V_400B),
                ('bCompressed'              ,'l'    ,V_400B),
                ('Unknown2'                 ,'1422B',V_400B),
                ('hExpectedPaddings'        ,'h'    ,V_430)
            )}

    @property
    def __h_elts(self):
        return self.__h_elt_versions[self.__version_bin]


class UnknownPaddingHeader(BiopacHeader):
    """
    I don't know what this is for, but it's 40-bytes long and right before some
    modern files
    """
    def __init__(self, file_revision, byte_order_char, **kwargs):
        self.file_revision = file_revision
        super().__init__(
            self.__h_elts, file_revision, byte_order_char, **kwargs)

    @property
    def effective_len_bytes(self):
        return self.data['lChannelLen']

    @property
    def __h_elts(self):
        elts = VersionedHeaderStructure(
            ('lChannelLen', 'l', V_ALL),
            ('Uknown', '36B', V_ALL)
        )
        return elts


class ChannelHeader(BiopacHeader):
    """
    The main Channel Header for an AcqKnowledge file. Note that this is known
    to be wrong for more modern files -- but for the purposes of this
    reader, we don't care. Enough fields are right that we can still find
    our way around.
    """

    def __init__(self, file_revision, byte_order_char, **kwargs):
        self.file_revision = file_revision
        super().__init__(self.__h_elts, file_revision, byte_order_char,
                         **kwargs)

    @property
    def __version_bin(self):
        bin = 'Unknown'
        if self.file_revision < V_400B:
            bin = 'PRE_4'
        else:
            bin = 'POST_4'
        return bin

    @property
    def effective_len_bytes(self):
        return self.data['lChanHeaderLen']

    @property
    def frequency_divider(self):
        """
        The number you divide the graph's samples_per_second by to get this
        channel's samples_per_second. If not defined (as in old file
        versions), this should be 1.
        """
        return self.data.get('nVarSampleDivider') or 1

    @property
    def raw_scale(self):
        return self.data['dAmplScale']

    @property
    def raw_offset(self):
        return self.data['dAmplOffset']

    @property
    def units(self):
        return self.data['szUnitsText'].decode(
            self.encoding, errors='ignore').strip('\0')

    @property
    def name(self):
        return self.data['szCommentText'].decode(
            self.encoding, errors='ignore').strip('\0')

    @property
    def point_count(self):
        return self.data['lBufLength']

    @property
    def order_num(self):
        return self.data['nChanOrder']

    @property
    def __h_elts(self):
        return self.__h_elt_versions[self.__version_bin]

    @property
    def __h_elt_versions(self):
        return {
            'PRE_4' : VersionedHeaderStructure(
                ('lChanHeaderLen'           ,'l'    ,V_20a),
                ('nNum'                     ,'h'    ,V_20a),
                ('szCommentText'            ,'40s'  ,V_20a),
                ('rgbColor'                 ,'4B'   ,V_20a),
                ('nDispChan'                ,'h'    ,V_20a),
                ('dVoltOffset'              ,'d'    ,V_20a),
                ('dVoltScale'               ,'d'    ,V_20a),
                ('szUnitsText'              ,'20s'  ,V_20a),
                ('lBufLength'               ,'l'    ,V_20a),
                ('dAmplScale'               ,'d'    ,V_20a),
                ('dAmplOffset'              ,'d'    ,V_20a),
                ('nChanOrder'               ,'h'    ,V_20a),
                ('nDispSize'                ,'h'    ,V_20a),
                ('plotMode'                 ,'h'    ,V_30r),
                ('vMid'                     ,'d'    ,V_30r),
                ('szDescription'            ,'128s' ,V_370),
                ('nVarSampleDivider'        ,'h'    ,V_370),
                ('vertPrecision'            ,'h'    ,V_373),
                ('activeSegmentColor'       ,'4b'   ,V_382),
                ('activeSegmentStyle'       ,'l'    ,V_382),
            ),
            'POST_4' : VersionedHeaderStructure(
                ('lChanHeaderLen'           ,'l'    ,V_20a),
                ('nNum'                     ,'h'    ,V_20a),
                ('szCommentText'            ,'40s'  ,V_20a),
                ('notColor'                 ,'4B'   ,V_20a),
                ('nDispChan'                ,'h'    ,V_20a),
                ('dVoltOffset'              ,'d'    ,V_20a),
                ('dVoltScale'               ,'d'    ,V_20a),
                ('szUnitsText'              ,'20s'  ,V_20a),
                ('lBufLength'               ,'l'    ,V_20a),
                ('dAmplScale'               ,'d'    ,V_20a),
                ('dAmplOffset'              ,'d'    ,V_20a),
                ('nChanOrder'               ,'h'    ,V_20a),
                ('nDispSize'                ,'h'    ,V_20a),
                ('unknown'                  ,'40s'  ,V_400B),
                ('nVarSampleDivider'        ,'h'    ,V_400B),
            )}


class ForeignHeader(BiopacHeader):
    """
    The Foreign Data header for AcqKnowledge files. I admit that I don't
    know what goes in this section.
    """

    def __init__(self, file_revision, byte_order_char, **kwargs):
        self.file_revision = file_revision
        super().__init__(self.__h_elts, file_revision, byte_order_char,
                         **kwargs)

    @property
    def __version_bin(self):
        bin = 'Unknown'
        if self.file_revision <= V_390:
            bin = "PRE_4"
        elif self.file_revision >= V_400B:
            bin = "POST_4"
        return bin

    @property
    def effective_len_bytes(self):
        if self.__version_bin == "PRE_4":
            return self.data['nLength']
        return self.data['lLength']

    @property
    def __h_elts(self):
        return self.__h_elt_versions[self.__version_bin]

    @property
    def __h_elt_versions(self):
        return {
            "PRE_4" : VersionedHeaderStructure(
                ('nLength'                  ,'h'    ,V_20a),
                ('nType'                    ,'h'    ,V_20a),
            ),
            "POST_4" : VersionedHeaderStructure(
                ('lLength'                  ,'l'    ,V_400B),
            )}


class ChannelDTypeHeader(BiopacHeader):
    def __init__(self, file_revision, byte_order_char, **kwargs):
        super().__init__(self.__h_elts, file_revision, byte_order_char,
                         **kwargs)

    @property
    def type_code(self):
        return self.data['nType']

    CODE_MAP = {
        0: 'f8',
        1: 'f8',
        2: 'i2'
    }

    @property
    def possibly_valid(self):
        # TODO: When and if we someday figure out the post-foreign-data header
        # mystery, remove this. Until then, we'll use this function to guess at
        # the location of the channel dtype header set.
        dtype_code = self.CODE_MAP.get(self.type_code, None)
        if dtype_code is None:
            return False
        type_size = int(dtype_code[-1])
        return type_size == self.sample_size

    @property
    def numpy_dtype(self):
        """
        A valid numpy dtype string: byte order and data type.

        Examples: '<i2', '>f8'
        """
        return self.byte_order_char + self.CODE_MAP[self.type_code]

    @property
    def sample_size(self):
        return self.data['nSize']

    @property
    def __h_elts(self):
        # This lets the standard effective_len_bytes work fine, I think.
        return VersionedHeaderStructure(
        ('nSize'                    ,'h'    ,V_20a),
        ('nType'                    ,'h'    ,V_20a),
        )


class V2JournalHeader(BiopacHeader):
    """
    Version 2-3 journal headers are trivial -- there's a four-byte tag that
    always contains 0x44332211, followed by a boolean "show" and then the
    length of the journal text.
    """
    EXPECTED_TAG_VALUE = (0x44, 0x33, 0x22, 0x11)
    EXPECTED_TAG_VALUE_HEX = "".join(f"{b:02X}" for b in EXPECTED_TAG_VALUE)

    def __init__(self, file_revision, byte_order_char, **kwargs):
        super().__init__(self.__h_elts, file_revision, byte_order_char,
                         **kwargs)

    @property
    def __h_elts(self):
        return VersionedHeaderStructure(
            ('tag',         '4B',   V_20a),
            ('hShow',       'h',    V_20a),
            ('lJournalLen', 'l',    V_20a)
        )

    @property
    def show(self):
        return self.data['hShow']

    @property
    def journal_len(self):
        return self.data['lJournalLen']
    
    @property
    def tag_value(self):
        return self.data['tag']
    
    @property
    def tag_value_hex(self):
        return "".join(f"{b:02X}" for b in self.data['tag'])
    
    def tag_value_matches_expected(self):
        return self.tag_value == self.EXPECTED_TAG_VALUE


class V4JournalLengthHeader(BiopacHeader):
    """
    In the case where there's no journal data, there's no full journal header.
    Instead, we just have a single long that tells us how much journal stuff
    (data + header) there is. Basically, if this value is less than the
    length of the V4JournalHeader, don't even try to read that header or
    journal data.

    The next stuff (if there is any) will be at self.offset + lJournalDataLen
    """

    def __init__(self, file_revision, byte_order_char, **kwargs):
        super().__init__(self.__h_elts, file_revision, byte_order_char,
                         **kwargs)

    @property
    def __h_elts(self):
        return VersionedHeaderStructure(
            ('lJournalDataLen', 'l', V_400B)
        )

    @property
    def journal_len(self):
        return self.data['lJournalDataLen']

    @property
    def data_end(self):
        return self.offset + self.journal_len



class V4JournalHeader(BiopacHeader):
    """
    In Version 4.1 and less, the journal is stored as plain text. From 4.2,
    it's stored as HTML. The start of the header tells the length of the
    entire journal section -- journal text and some preamble; the compression
    headers (if compressed) follow at self.offset + lFullLength.

    The length of the actual journal text is contained 266 bytes into the
    preamble for 4.1 files; in 4.2 files it's a long immediately before the
    journal data.

    The journal data starts at 560 bytes after the start of this header in
    4.1 files; in 4.2 it's at 594 bytes, and 4.4 it's at 598 bytes.
    """

    def __init__(self, file_revision, byte_order_char, **kwargs):
        super().__init__(self.__h_elts, file_revision, byte_order_char,
                         **kwargs)

    @property
    def __h_elts(self):
        return VersionedHeaderStructure(
            ('bUnknown1', '262b', V_400B),
            ('lEarlyJournalLen', 'l', V_400B),
            ('bUnknown2', '290b', V_400B),
            ('bUnknown3', '26b', V_420),
            ('bUnknown4', '4b', V_440),
            ('lLateJournalLenMinusOne', 'l', V_420),
            ('lLateJournalLen', 'l', V_420)
        )

    @property
    def journal_len(self):
        if self.file_revision < V_420:
            return self.data['lEarlyJournalLen']
        return self.data['lLateJournalLen']


class MainCompressionHeader(BiopacHeader):
    # In compressed files, the markers are stored where the data would be in
    # uncompressed files. There's also some padding, and I don't know
    # what's in it.
    def __init__(self, file_revision, byte_order_char, **kwargs):
        self.file_revision = file_revision
        super().__init__(self.__h_elts, file_revision, byte_order_char,
                         **kwargs)

    @property
    def effective_len_bytes(self):
        return self.__effective_len_byte_versions[self.__version_bin]()

    @property
    def __effective_len_byte_versions(self):
        # Determined through experimentation, may not be correct for some
        # revisions. Or files, for that matter. We'll see.
        return {
            'PRE_4': lambda: (
                self.struct_dict.len_bytes + self.data['lTextLen']),
            'POST_4': lambda: (  # This is incorrect for data with journals
                self.struct_dict.len_bytes +
                self.data['lStrLen1'] +
                self.data['lStrLen2']
                )
        }

    @property
    def __version_bin(self):
        bin = "Unknown"
        if self.file_revision <= V_400B:
            bin = "PRE_4"
        else:
            bin = "POST_4"
        return bin

    @property
    def __h_elts(self):
        return self.__h_elts_versions[self.__version_bin]

    @property
    def __h_elts_versions(self):
        return {
            'PRE_4': VersionedHeaderStructure(
                ('Unknown', '34B', V_20a),
                ('lTextLen', 'l', V_20a)
            ),
            'POST_4': VersionedHeaderStructure(
                ('Unknown1', '24B', V_400B),  # Should probably be 24.
                ('lStrLen1', 'l', V_400B),
                ('lStrLen2', 'l', V_400B),
                ('Unknown2', '20B', V_400B),
                ('Unknown3', '6B', V_420),
            )
        }


class ChannelCompressionHeader(BiopacHeader):
    def __init__(self, file_revision, byte_order_char, **kwargs):
        super().__init__(self.__h_elts, file_revision, byte_order_char,
                         **kwargs)

    @property
    def effective_len_bytes(self):
        """
        Return the length of the header UP TO THE NEXT HEADER, skipping the
        compressed data. Use header_only_len_bytes for only the header length.
        """
        return self.header_only_len_bytes + self.data['lCompressedLen']

    @property
    def header_only_len_bytes(self):
        """
        A truly variable-length header. Oddly, the channel description and
        units are included in the header. This class doesn't have a way to
        get those out, but we'll just use them from the earlier channel
        header.
        Immediately after this header, the compressed data starts -- it will
        be after the units text, and starts with 'x'.
        """
        return (
            self.struct_dict.len_bytes + self.data['lChannelLabelLen'] +
            self.data['lUnitLabelLen'])

    @property
    def compressed_data_offset(self):
        """
        The offset to the compressed data.
        Note that this won't be valid until self#unpack_from_file() is run.
        """
        return self.offset + self.header_only_len_bytes

    @property
    def compressed_data_len(self):
        return self.data['lCompressedLen']

    @property
    def __h_elts(self):
        return VersionedHeaderStructure(
        ('Unknown'                  ,'44B'  ,V_381),
        ('lChannelLabelLen'         ,'l'    ,V_381),
        ('lUnitLabelLen'            ,'l'    ,V_381),
        ('lUncompressedLen'         ,'l'    ,V_381),
        ('lCompressedLen'           ,'l'    ,V_381),
        )


class V2MarkerHeader(BiopacHeader):
    """
    Marker structure for files in Version 3, very likely down to version 2.
    """

    def __init__(self, file_revision, byte_order_char, **kwargs):
        super().__init__(self.__h_elts, file_revision, byte_order_char,
                         **kwargs)

    # NOTE: lLength does NOT include this header length -- only the length
    # of all the marker items. This is different than in the v4 header.
    @property
    def __h_elts(self):
        return VersionedHeaderStructure(
        ('lLength'              ,'l'    ,V_20a),
        ('lMarkers'             ,'l'    ,V_20a),
        )

    @property
    def marker_count(self):
        return self.data['lMarkers']

class V2MarkerMetadataPreHeader(BiopacHeader):
    """
    Tells us how many marker metadata headers there are. It's probably the same
    as the marker count, but it's stored separately.
    """

    def __init__(self, file_revision, byte_order_char, **kwargs):
        super().__init__(self.__h_elts, file_revision, byte_order_char,
                         **kwargs)

    @property
    def __h_elts(self):
        return VersionedHeaderStructure(
        ('tag'           ,'4B'   ,V_20a),
        ('lItemCount'    ,'l'    ,V_20a),
        ('sUnknown'      ,'76s'  ,V_20a),
        )
    
    @property
    def item_count(self):
        return self.data['lItemCount']
    
    @property
    def tag_value(self):
        return self.data['tag']
    
    
class V2MarkerMetadataHeader(BiopacHeader):
    """
    Marker metadata for files in Version 2.
    """

    def __init__(self, file_revision, byte_order_char, **kwargs):
        super().__init__(self.__h_elts, file_revision, byte_order_char,
                         **kwargs)  

    @property
    def __h_elts(self):
        return VersionedHeaderStructure(
        ('lUnknown1'            ,'l'    ,V_20a),
        ('lMarkerNumber'        ,'l'    ,V_20a),
        ('bUnknown2'            ,'12B'  ,V_20a),
        ('rgbaColor'            ,'4B'   ,V_20a),
        ('hMarkerTag'           ,'h'    ,V_20a),
        ('hMarkerTypeId'        ,'h'    ,V_20a),
        )
    
    @property
    def marker_number(self):
        return self.data['lMarkerNumber']
    
    @property
    def rgba_color(self):
        return self.data['rgbaColor']
    
    @property
    def marker_tag(self):
        return self.data['hMarkerTag']    
    
    @property
    def marker_index(self):
        return self.marker_number - 1
    

class V4MarkerHeader(BiopacHeader):
    """
    Marker structure for files from Version 4 onwards
    """

    def __init__(self, file_revision, byte_order_char, **kwargs):
        super().__init__(self.__h_elts, file_revision, byte_order_char,
                         **kwargs)

    # NOTE: lLength INCLUDES this header length -- the markers end at
    # marker_start_offset + lLength. This is different than in the v2 header.
    @property
    def __h_elts(self):
        return VersionedHeaderStructure(
        ('lLength'              ,'l'    ,V_400B),
        ('lMarkersExtra'        ,'l'    ,V_400B),
        ('lMarkers'             ,'l'    ,V_400B),
        ('Unknown'              ,'6B'   ,V_400B),
        ('szDefl'               ,'5s'   ,V_400B),
        ('Unknown2'             ,'h'    ,V_400B),
        ('Unknown3'             ,'8B'   ,V_42x),
        ('Unknown4'             ,'8B'   ,V_440)
        )

    # I'm not quite sure about these two marker count headers; they seem to
    # both be wrong.
    @property
    def marker_count(self):
        return self.data['lMarkersExtra'] - 1


class V2MarkerItemHeader(BiopacHeader):
    """
    Marker Items for files in Version 3, very likely down to version 2.
    """

    def __init__(self, file_revision, byte_order_char, **kwargs):
        super().__init__(self.__h_elts, file_revision, byte_order_char,
                         **kwargs)

    @property
    def __h_elts(self):
        return VersionedHeaderStructure(
        ('lSample'              ,'l'    ,V_20a),
        ('fSelected'            ,'h'    ,V_35x),
        ('fTextLocked'          ,'h'    ,V_20a),
        ('fPositionLocked'      ,'h'    ,V_20a),
        ('nTextLength'          ,'h'    ,V_20a),
        )

    # Note: The spec says nTextLength includes the trailing null, but it
    # seems to not...?
    # It seems that it does include it in V_303 so haha
    @property
    def text_length(self):
        if self.file_revision < V_35x:
            return self.data['nTextLength']
        else:
            return self.data['nTextLength'] + 1

    @property
    def sample_index(self):
        return self.data['lSample']

    @property
    def channel_number(self):
        """ None means it's a global marker """
        return None

    @property
    def date_created_ms(self):
        """ markers don't get "Date created" until v. 4.4.0 """
        return None

    @property
    def type_code(self):
        """ These markers don't get type_codes, but it's OK """
        return None



class V4MarkerItemHeader(BiopacHeader):
    """
    Marker Items for files in Version 4 ownards.
    """

    def __init__(self, file_revision, byte_order_char, **kwargs):
        super().__init__(self.__h_elts, file_revision, byte_order_char,
                         **kwargs)

    @property
    def __h_elts(self):
        return VersionedHeaderStructure(
        ('lSample'              ,'l'    ,V_400B),
        ('Unknown'              ,'4B'   ,V_400B),
        ('nChannel'             ,'h'    ,V_400B),
        ('sMarkerStyle'         ,'4s'   ,V_400B),
        ('llDateCreated'        ,'Q'    ,V_440),
        ('Unknown3'             ,'8B'   ,V_42x),
        ('nTextLength'          ,'h'    ,V_400B),
        )

    # Unlike in older versions, nTextLength does include the trailing null.
    @property
    def text_length(self):
        return self.data['nTextLength']

    @property
    def sample_index(self):
        return self.data['lSample']

    @property
    def channel_number(self):
        """ None means it's a global marker """
        chan = self.data['nChannel']
        if chan == -1:
            chan = None
        return chan

    @property
    def date_created_ms(self):
        """ Date when marker was created (in ms since 1970-01-01) """
        if self.file_revision < V_440:
            return None
        else:
            return self.data['llDateCreated']

    @property
    def type_code(self):
        return self.data['sMarkerStyle'].decode(self.encoding, errors='ignore')
