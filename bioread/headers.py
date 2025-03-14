# coding: utf8
# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2025 Board of Regents of the University of Wisconsin System
#
# Written by Nate Vack <njvack@wisc.edu> with research from John Ollinger
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison, and later the Center for Healthy Minds

# flake8: noqa: E203, E231, E741, E405, E403

from bioread.file_revisions import *

import ctypes


class HeaderDefinitionError(TypeError):
    """
    Exception raised when a header definition is invalid.
    """

    pass


class HeaderMeta(type):
    """
    Metaclass for the Header class. We'll use this to register subclasses
    that are valid for a given revision, so you can do something like:
    main_compression_header = ChannelHeader.for_revision(file_revision, ...)
    and you'll get either a ChannelHeaderPre4 or ChannelHeaderPost4 depenting
    on file_revision.

    Note that to register a class as a subclass, you need to both
    be a subclass of Header and have an is_valid_for_revision() method.
    """

    def __new__(mcs, name, bases, class_dict):
        # Initialize _versioned_subclasses for this class if it doesn't exist
        cls = super().__new__(mcs, name, bases, class_dict)
        # Don't register Header itself, it's an abstract base class
        if name == "Header":
            return cls
        if "_versioned_subclasses" not in class_dict:
            cls._versioned_subclasses = []
        # Only register if this class has its own is_valid_for_revision method
        if "is_valid_for_revision" in class_dict:
            # Register on the direct parent class
            superclass = bases[0]
            if not hasattr(superclass, "_versioned_subclasses"):
                superclass._versioned_subclasses = []
            superclass._versioned_subclasses.append(cls)
        return cls


def revision_range(min_revision=None, max_revision=None):
    """
    Create a class method that checks if a revision is within the given range.
    """

    def is_valid_for_revision(cls, file_revision):
        if min_revision is not None and file_revision < min_revision:
            return False
        if max_revision is not None and file_revision >= max_revision:
            return False
        return True

    return classmethod(is_valid_for_revision)


def pre_4():
    return revision_range(max_revision=V_400B)


def post_4():
    return revision_range(min_revision=V_400B)


class Header(metaclass=HeaderMeta):
    """
    The base class for all the Biopac headers we're going to use. Basically,
    this contains tooling to build a structure with appropriate endianness,
    handle the thing where fields appear in different versions of the AcqKnowledge,
    and then has the ability to unpack the data from a file.
    """

    BASES = {">": ctypes.BigEndianStructure, "<": ctypes.LittleEndianStructure}

    def __init__(
        self, file_revision, byte_order_char, encoding="utf-8", _from_factory=False
    ):
        """
        Private constructor. Use for_revision() to create a header object.
        If you really need to use this method, pass _from_factory=True.
        """
        if not _from_factory:
            raise ValueError("Headers should be created using for_revision() method.")
        self.file_revision = file_revision
        self.byte_order_char = byte_order_char
        self.encoding = encoding
        self.offset = None
        self.raw_data = None
        self._struct = None

    # @classmethod
    # def _register_subclass(cls, subclass):
    #     print(f"Registering subclass {subclass.__name__} on {cls.__name__}")
    #     cls._versioned_subclasses.append(subclass)

    @classmethod
    def for_revision(cls, file_revision, byte_order_char, encoding="utf-8"):
        """
        Factory method for creating a header object for a given revision.
        Checks _versioned_subclasses for a subclass that returns True for
        is_valid_for_revision(), and if it finds one, instantiates that.
        Otherwise, it returns an instance of the base class.
        """
        constructor_args = [file_revision, byte_order_char, encoding, True]
        matching_subclasses = [
            subclass
            for subclass in cls._versioned_subclasses
            if subclass.is_valid_for_revision(file_revision)
        ]
        if len(matching_subclasses) == 0:
            return cls(*constructor_args)

        if len(matching_subclasses) > 1:
            # This should never happen, but it's good to have this check
            raise HeaderDefinitionError(
                f"Multiple subclasses of {cls.__name__} are valid for revision {file_revision}: {matching_subclasses}"
            )

        return matching_subclasses[0](*constructor_args)

    def unpack_from_file(self, data_file, offset):
        """Unpack header data from a file at the given offset"""
        self.offset = offset
        data_file.seek(offset)
        # Read enough bytes for this structure
        self.unpack_from_bytes(data_file.read(self.struct_length))

    def unpack_from_bytes(self, bytes_data):
        """Unpack header data from a string/bytes. Useful for testing / debugging."""
        self.raw_data = bytes_data
        self._struct = self._struct_class.from_buffer_copy(bytes_data)

    @property
    def _struct_class(self):
        """
        Get the appropriate structure class with correct byte order
        """
        if not hasattr(self, "__cached_struct_class"):
            base_class = self.BASES[self.byte_order_char]
            fields = self._fields_for_current_version()

            class StructClass(base_class):
                _pack_ = 1
                _fields_ = fields

            self.__cached_struct_class = StructClass
        return self.__cached_struct_class

    @property
    def struct_length(self):
        """
        Length of the structure itself. You should probably not override this!
        """
        return ctypes.sizeof(self._struct_class)

    @property
    def effective_len_bytes(self):
        """
        The number of bytes it'll take to get to the next header in the file.
        This may use the information from the header -- generally in cases where
        the header contains things like variable length strings.
        """
        return self.struct_length

    def _fields_for_current_version(self):
        """
        Return a list of fields that should be included in the structure.
        """
        return [
            (field[0], field[1])
            for field in self.__class__._versioned_fields
            if field[2] <= self.file_revision
        ]

    @property
    def data(self):
        """
        Return a dictionary of the fields in the structure.
        """
        if not hasattr(self, "__data"):
            self.__data = {}
            for field, _field_type in self._struct._fields_:
                if isinstance(getattr(self._struct, field), ctypes.Array):
                    self.__data[field] = list(getattr(self._struct, field))
                else:
                    self.__data[field] = getattr(self._struct, field)
        return self.__data


class GraphHeader(Header):
    """
    Base class for Graph Headers to yield a common interface.
    """

    @property
    def effective_len_bytes(self):
        return self._struct.lExtItemHeaderLen

    @property
    def channel_count(self):
        return self._struct.nChannels

    @property
    def sample_time(self):
        return self._struct.dSampleTime

    @property
    def compressed(self):
        if not hasattr(self._struct, "bCompressed"):
            return False
        return self._struct.bCompressed != 0

    @property
    def data_format(self):
        return "compressed" if self.compressed else "uncompressed"

    @property
    def expected_padding_headers(self):
        if self.file_revision >= V_430:
            return self._struct.hExpectedPaddings
        return 0


class GraphHeaderPre4(GraphHeader):
    """
    Graph Header for files with revision less than 4.
    """

    is_valid_for_revision = pre_4()

    _versioned_fields = [
        ("nItemHeaderLen", ctypes.c_int16, V_ALL),
        ("lVersion", ctypes.c_uint32, V_ALL),
        ("lExtItemHeaderLen", ctypes.c_int32, V_20a),
        ("nChannels", ctypes.c_int16, V_20a),
        ("nHorizAxisType", ctypes.c_int16, V_20a),
        ("nCurChannel", ctypes.c_int16, V_20a),
        ("dSampleTime", ctypes.c_double, V_20a),
        ("dTimeOffset", ctypes.c_double, V_20a),
        ("dTimeScale", ctypes.c_double, V_20a),
        ("dTimeCursor1", ctypes.c_double, V_20a),
        ("dTimeCursor2", ctypes.c_double, V_20a),
        ("rcWindow", ctypes.c_int16 * 4, V_20a),
        ("nMeasurement", ctypes.c_int16 * 6, V_20a),
        ("fHilite", ctypes.c_int16, V_20a),
        ("dFirstTimeOffset", ctypes.c_double, V_20a),
        ("nRescale", ctypes.c_int16, V_20a),
        ("szHorizUnits1", ctypes.c_char * 40, V_20a),
        ("szHorizUnits2", ctypes.c_char * 10, V_20a),
        ("nInMemory", ctypes.c_int16, V_20a),
        ("fGrid", ctypes.c_int16, V_20a),
        ("fMarkers", ctypes.c_int16, V_20a),
        ("nPlotDraft", ctypes.c_int16, V_20a),
        ("nDispMode", ctypes.c_int16, V_20a),
        ("rRReserved", ctypes.c_int16, V_20a),
        ("BShowToolBar", ctypes.c_int16, V_30r),
        ("BShowChannelButtons", ctypes.c_int16, V_30r),
        ("BShowMeasurements", ctypes.c_int16, V_30r),
        ("BShowMarkers", ctypes.c_int16, V_30r),
        ("BShowJournal", ctypes.c_int16, V_30r),
        ("CurXChannel", ctypes.c_int16, V_30r),
        ("MmtPrecision", ctypes.c_int16, V_30r),
        ("NMeasurementRows", ctypes.c_int16, V_303),
        ("mmt40", ctypes.c_int16 * 40, V_303),
        ("mmtChan40", ctypes.c_int16 * 40, V_303),
        ("MmtCalcOpnd1", ctypes.c_int16 * 40, V_35x),
        ("MmtCalcOpnd2", ctypes.c_int16 * 40, V_35x),
        ("MmtCalcOp", ctypes.c_int16 * 40, V_35x),
        ("MmtCalcConstant", ctypes.c_double * 40, V_35x),
        ("bNewGridWithMinor", ctypes.c_int32, V_370),
        ("colorMajorGrid", ctypes.c_ubyte * 4, V_370),
        ("colorMinorGrid", ctypes.c_ubyte * 4, V_370),
        ("wMajorGridStyle", ctypes.c_uint16, V_370),
        ("wMinorGridStyle", ctypes.c_uint16, V_370),
        ("wMajorGridWidth", ctypes.c_uint16, V_370),
        ("wMinorGridWidth", ctypes.c_uint16, V_370),
        ("bFixedUnitsDiv", ctypes.c_uint32, V_370),
        ("bMid_Range_Show", ctypes.c_uint32, V_370),
        ("dStart_Middle_Point", ctypes.c_double, V_370),
        ("dOffset_Point", ctypes.c_double * 60, V_370),
        ("hGrid", ctypes.c_double, V_370),
        ("vGrid", ctypes.c_double * 60, V_370),
        ("bEnableWaveTools", ctypes.c_int32, V_370),
        ("hozizPrecision", ctypes.c_int16, V_373),
        ("Reserved", ctypes.c_byte * 20, V_381),
        ("bOverlapMode", ctypes.c_int32, V_381),
        ("bShowHardware", ctypes.c_int32, V_381),
        ("bXAutoPlot", ctypes.c_int32, V_381),
        ("bXAutoScroll", ctypes.c_int32, V_381),
        ("bStartButtonVisible", ctypes.c_int32, V_381),
        ("bCompressed", ctypes.c_int32, V_381),
        ("bAlwaysStartButtonVisible", ctypes.c_int32, V_381),
        ("pathVideo", ctypes.c_char * 260, V_382),
        ("optSyncDelay", ctypes.c_int32, V_382),
        ("syncDelay", ctypes.c_double, V_382),
        ("bHRP_PasteMeasurements", ctypes.c_int32, V_382),
        ("graphType", ctypes.c_int32, V_390),
        ("mmtCalcExpr", ctypes.c_char * 10240, V_390),
        ("mmtMomentOrder", ctypes.c_int32 * 40, V_390),
        ("mmtTimeDelay", ctypes.c_int32 * 40, V_390),
        ("mmtEmbedDim", ctypes.c_int32 * 40, V_390),
        ("mmtMIDelay", ctypes.c_int32 * 40, V_390),
    ]


class GraphHeaderPost4(GraphHeader):
    """
    Graph Header for files with revision 4 and above.
    """

    is_valid_for_revision = post_4()

    _versioned_fields = [
        ("nItemHeaderLen", ctypes.c_int16, V_ALL),
        ("lVersion", ctypes.c_int32, V_ALL),
        ("lExtItemHeaderLen", ctypes.c_int32, V_20a),
        ("nChannels", ctypes.c_int16, V_20a),
        ("nHorizAxisType", ctypes.c_int16, V_20a),
        ("nCurChannel", ctypes.c_int16, V_20a),
        ("dSampleTime", ctypes.c_double, V_20a),
        ("dTimeOffset", ctypes.c_double, V_20a),
        ("dTimeScale", ctypes.c_double, V_20a),
        ("dTimeCursor1", ctypes.c_double, V_20a),
        ("dTimeCursor2", ctypes.c_double, V_20a),
        ("rcWindow", ctypes.c_int16 * 4, V_20a),
        ("nMeasurement", ctypes.c_int16 * 6, V_20a),
        ("fHilite", ctypes.c_int16, V_20a),
        ("dFirstTimeOffset", ctypes.c_double, V_20a),
        ("nRescale", ctypes.c_int16, V_20a),
        ("szHorizUnits1", ctypes.c_char * 40, V_20a),
        ("szHorizUnits2", ctypes.c_char * 10, V_20a),
        ("nInMemory", ctypes.c_int16, V_20a),
        ("fGrid", ctypes.c_int16, V_20a),
        ("fMarkers", ctypes.c_int16, V_20a),
        ("nPlotDraft", ctypes.c_int16, V_20a),
        ("nDispMode", ctypes.c_int16, V_20a),
        ("rRReserved", ctypes.c_int16, V_20a),
        ("Unknown", ctypes.c_byte * 822, V_400B),
        ("bCompressed", ctypes.c_int32, V_400B),
        ("Unknown2", ctypes.c_byte * 1422, V_400B),
        ("hExpectedPaddings", ctypes.c_int16, V_430),
    ]


class UnknownPaddingHeader(Header):
    """
    I don't know what this is for, but it's 40-bytes long and right before some
    modern files
    """

    _versioned_fields = [
        ("lChannelLen", ctypes.c_int32, V_ALL),
        ("Uknown", ctypes.c_byte * 36, V_ALL),
    ]

    @property
    def effective_len_bytes(self):
        return self._struct.lChannelLen


class ChannelHeader(Header):
    """
    Base class for Channel Headers to yield a common interface.
    """

    @property
    def effective_len_bytes(self):
        return self._struct.lChanHeaderLen

    @property
    def frequency_divider(self):
        if hasattr(self._struct, "nVarSampleDivider"):
            return self._struct.nVarSampleDivider or 1
        else:
            return 1

    @property
    def raw_scale(self):
        return self._struct.dAmplScale

    @property
    def raw_offset(self):
        return self._struct.dAmplOffset

    @property
    def units(self):
        return self._struct.szUnitsText.decode(self.encoding, errors="ignore").strip(
            "\0"
        )

    @property
    def name(self):
        return self._struct.szCommentText.decode(self.encoding, errors="ignore").strip(
            "\0"
        )

    @property
    def point_count(self):
        return self._struct.lBufLength

    @property
    def order_num(self):
        return self._struct.nChanOrder


class ChannelHeaderPre4(ChannelHeader):
    """
    Channel Header for files with revision less than 4.
    """

    is_valid_for_revision = pre_4()

    _versioned_fields = [
        ("lChanHeaderLen", ctypes.c_int32, V_20a),
        ("nNum", ctypes.c_int16, V_20a),
        ("szCommentText", ctypes.c_char * 40, V_20a),
        ("rgbColor", ctypes.c_ubyte * 4, V_20a),
        ("nDispChan", ctypes.c_int16, V_20a),
        ("dVoltOffset", ctypes.c_double, V_20a),
        ("dVoltScale", ctypes.c_double, V_20a),
        ("szUnitsText", ctypes.c_char * 20, V_20a),
        ("lBufLength", ctypes.c_int32, V_20a),
        ("dAmplScale", ctypes.c_double, V_20a),
        ("dAmplOffset", ctypes.c_double, V_20a),
        ("nChanOrder", ctypes.c_int16, V_20a),
        ("nDispSize", ctypes.c_int16, V_20a),
        ("plotMode", ctypes.c_int16, V_30r),
        ("vMid", ctypes.c_double, V_30r),
        ("szDescription", ctypes.c_char * 128, V_370),
        ("nVarSampleDivider", ctypes.c_int16, V_370),
        ("vertPrecision", ctypes.c_int16, V_373),
        ("activeSegmentColor", ctypes.c_ubyte * 4, V_382),
        ("activeSegmentStyle", ctypes.c_int32, V_382),
    ]


class ChannelHeaderPost4(ChannelHeader):
    """
    Channel Header for files with revision 4 and above.
    """

    is_valid_for_revision = post_4()

    _versioned_fields = [
        ("lChanHeaderLen", ctypes.c_int32, V_20a),
        ("nNum", ctypes.c_int16, V_20a),
        ("szCommentText", ctypes.c_char * 40, V_20a),
        ("notColor", ctypes.c_ubyte * 4, V_20a),
        ("nDispChan", ctypes.c_int16, V_20a),
        ("dVoltOffset", ctypes.c_double, V_20a),
        ("dVoltScale", ctypes.c_double, V_20a),
        ("szUnitsText", ctypes.c_char * 20, V_20a),
        ("lBufLength", ctypes.c_int32, V_20a),
        ("dAmplScale", ctypes.c_double, V_20a),
        ("dAmplOffset", ctypes.c_double, V_20a),
        ("nChanOrder", ctypes.c_int16, V_20a),
        ("nDispSize", ctypes.c_int16, V_20a),
        ("unknown", ctypes.c_ubyte * 40, V_400B),
        ("nVarSampleDivider", ctypes.c_int16, V_400B),
    ]


class ForeignHeader(Header):
    """
    Abstract base class for foreign headers.
    I'm genuinely not sure what the foreign data is for.
    """

    pass


class ForeignHeaderPre4(ForeignHeader):
    is_valid_for_revision = pre_4()

    _versioned_fields = [
        ("nLength", ctypes.c_int16, V_20a),
        ("nType", ctypes.c_int16, V_20a),
    ]

    @property
    def effective_len_bytes(self):
        return self._struct.nLength


class ForeignHeaderPost4(ForeignHeader):
    is_valid_for_revision = post_4()

    _versioned_fields = [("lLength", ctypes.c_int32, V_400B)]

    @property
    def effective_len_bytes(self):
        return self._struct.lLength


class ChannelDTypeHeader(Header):
    """
    Channel data types. This comes immediately after the foreign data, but occasionally
    there seems to be extra data -- so we use possibly_valid and scan forward until
    we find valid channel data types, like some kind of animal.
    """

    _versioned_fields = [
        ("nSize", ctypes.c_int16, V_20a),
        ("nType", ctypes.c_int16, V_20a),
    ]

    @property
    def type_code(self):
        return self._struct.nType

    CODE_MAP = {0: "f8", 1: "f8", 2: "i2"}

    @property
    def possibly_valid(self):
        dtype_code = self.CODE_MAP.get(self.type_code, None)
        if dtype_code is None:
            return False
        type_size = int(dtype_code[-1])
        return type_size == self.sample_size

    def numpy_dtype(self, data_is_compressed):
        # Compressed data is always little-endian
        byte_order_char = "<" if data_is_compressed else self.byte_order_char
        return byte_order_char + self.CODE_MAP[self.type_code]

    @property
    def sample_size(self):
        return self._struct.nSize


class JournalHeader(Header):
    """
    Abstract base class for journal headers.
    """

    pass


class JournalHeaderPre4(JournalHeader):
    """
    Version 2-3 journal headers are trivial -- there's a four-byte tag that
    always contains 0x44332211, followed by a boolean "show" and then the
    length of the journal text.
    """

    is_valid_for_revision = pre_4()

    EXPECTED_TAG_VALUE = (0x44, 0x33, 0x22, 0x11)
    EXPECTED_TAG_VALUE_HEX = "".join(f"{b:02X}" for b in EXPECTED_TAG_VALUE)

    _versioned_fields = [
        ("tag", ctypes.c_ubyte * 4, V_20a),
        ("hShow", ctypes.c_int16, V_20a),
        ("lJournalLen", ctypes.c_int32, V_20a),
    ]

    @property
    def show(self):
        return self._struct.hShow

    @property
    def journal_len(self):
        return self._struct.lJournalLen

    @property
    def tag_value(self):
        return tuple(self._struct.tag)

    @property
    def tag_value_hex(self):
        return "".join(f"{b:02X}" for b in self._struct.tag)

    def tag_value_matches_expected(self):
        return self.tag_value == self.EXPECTED_TAG_VALUE


class JournalHeaderPost4(JournalHeader):
    """
    In Version 4.1 and less, the journal is stored as plain text. From 4.2,
    it's stored as HTML. The start of the header tells the length of the
    entire journal section -- journal text and some preamble; the compression
    headers (if compressed) follow at self.offset + lFullLength.
    """

    is_valid_for_revision = post_4()

    _versioned_fields = [
        ("bUnknown1", ctypes.c_byte * 262, V_400B),
        ("lEarlyJournalLen", ctypes.c_int32, V_400B),
        ("bUnknown2", ctypes.c_byte * 290, V_400B),
        ("bUnknown3", ctypes.c_byte * 26, V_420),
        ("bUnknown4", ctypes.c_byte * 4, V_440),
        ("lLateJournalLenMinusOne", ctypes.c_int32, V_420),
        ("lLateJournalLen", ctypes.c_int32, V_420),
    ]

    @property
    def journal_len(self):
        if self.file_revision < V_420:
            return self._struct.lEarlyJournalLen
        return self._struct.lLateJournalLen


class JournalLengthHeader(Header):
    """
    In the case where there's no journal data, there's no full journal header.
    Instead, we just have a single long that tells us how much journal stuff
    (data + header) there is. Basically, if this value is less than the
    length of the V4JournalHeader, don't even try to read that header or
    journal data.

    The next stuff (if there is any) will be at self.offset + lJournalDataLen

    This is only present in version 4.
    """

    # Define the structure fields
    _versioned_fields = [("lJournalDataLen", ctypes.c_int32, V_400B)]

    @property
    def journal_len(self):
        """Get the journal length directly from the structure"""
        return self._struct.lJournalDataLen

    @property
    def data_end(self):
        """Calculate the end position of the journal data"""
        return self.offset + self.journal_len


class MainCompressionHeader(Header):
    """
    Abstract base class for main compression headers.
    """

    pass


class MainCompressionHeaderPre4(MainCompressionHeader):
    is_valid_for_revision = pre_4()

    _versioned_fields = [
        ("Unknown", ctypes.c_byte * 34, V_20a),
        ("lTextLen", ctypes.c_int32, V_20a),
    ]

    @property
    def effective_len_bytes(self):
        return self.struct_length + self._struct.lTextLen


class MainCompressionHeaderPost4(MainCompressionHeader):
    is_valid_for_revision = post_4()

    _versioned_fields = [
        ("Unknown1", ctypes.c_byte * 24, V_400B),
        ("lStrLen1", ctypes.c_int32, V_400B),
        ("lStrLen2", ctypes.c_int32, V_400B),
        ("Unknown2", ctypes.c_byte * 20, V_400B),
        ("Unknown3", ctypes.c_byte * 6, V_420),
    ]

    @property
    def effective_len_bytes(self):
        return self.struct_length + self._struct.lStrLen1 + self._struct.lStrLen2


class ChannelCompressionHeader(Header):
    """
    Represents the channel compression header in an AcqKnowledge file.
    """

    _versioned_fields = [
        ("Unknown", ctypes.c_byte * 44, V_381),
        ("lChannelLabelLen", ctypes.c_int32, V_381),
        ("lUnitLabelLen", ctypes.c_int32, V_381),
        ("lUncompressedLen", ctypes.c_int32, V_381),
        ("lCompressedLen", ctypes.c_int32, V_381),
    ]

    @property
    def effective_len_bytes(self):
        """
        Return the length of the header UP TO THE NEXT HEADER, skipping the
        compressed data. Use header_only_len_bytes for only the header length.
        """
        return self.header_only_len_bytes + self._struct.lCompressedLen

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
            self.struct_length
            + self._struct.lChannelLabelLen
            + self._struct.lUnitLabelLen
        )

    @property
    def compressed_data_offset(self):
        """
        The offset to the compressed data.
        Note that this won't be valid until self#unpack_from_file() is run.
        """
        return self.offset + self.header_only_len_bytes

    @property
    def compressed_data_len(self):
        return self._struct.lCompressedLen

    @property
    def uncompressed_data_len(self):
        return self._struct.lUncompressedLen


class MarkerHeader(Header):
    """
    Abstract base class for marker headers.
    """

    pass


class MarkerHeaderPre4(MarkerHeader):
    """
    Marker structure for files in Version 3, very likely down to version 2.
    """

    is_valid_for_revision = pre_4()

    _versioned_fields = [
        ("lLength", ctypes.c_int32, V_20a),
        ("lMarkers", ctypes.c_int32, V_20a),
    ]

    @property
    def marker_count(self):
        return self._struct.lMarkers


class MarkerHeaderPost4(MarkerHeader):
    """
    Marker structure for files from Version 4 onwards
    """

    is_valid_for_revision = post_4()

    _versioned_fields = [
        ("lLength", ctypes.c_int32, V_400B),
        ("lMarkersExtra", ctypes.c_int32, V_400B),
        ("lMarkers", ctypes.c_int32, V_400B),
        ("Unknown", ctypes.c_byte * 6, V_400B),
        ("szDefl", ctypes.c_char * 5, V_400B),
        ("Unknown2", ctypes.c_int16, V_400B),
        ("Unknown3", ctypes.c_byte * 8, V_42x),
        ("Unknown4", ctypes.c_byte * 8, V_440),
    ]

    @property
    def marker_count(self):
        return self._struct.lMarkersExtra - 1


class MarkerPreItemMetadataHeaderPre4(Header):
    """
    I'm not sure what the data here mean -- there's an item count field, but
    the number from MarkerHeaderPre4 seems to be the correct one.
    """

    _versioned_fields = [
        ("tag", ctypes.c_byte * 4, V_20a),
        ("lItemCount", ctypes.c_int32, V_20a),
        ("sUnknown", ctypes.c_ubyte * 76, V_20a),
    ]

    @property
    def tag_value(self):
        # Convert the c_byte_Array_4 to a tuple
        return tuple(self._struct.tag)


class MarkerItemMetadataHeader(Header):
    """
    Marker metadata for files in Version 2-3.
    This header is not present in Version 4.
    """

    _versioned_fields = [
        ("lUnknown1", ctypes.c_int32, V_20a),
        ("lMarkerNumber", ctypes.c_int32, V_20a),
        ("bUnknown2", ctypes.c_byte * 12, V_20a),
        ("rgbaColor", ctypes.c_byte * 4, V_20a),
        ("hMarkerTag", ctypes.c_int16, V_20a),
        ("hMarkerTypeId", ctypes.c_int16, V_20a),
    ]

    @property
    def marker_number(self):
        return self._struct.lMarkerNumber

    @property
    def rgba_color(self):
        return tuple(self._struct.rgbaColor)

    @property
    def marker_tag(self):
        return self._struct.hMarkerTag

    @property
    def marker_index(self):
        return self.marker_number - 1


class MarkerItemHeader(Header):
    """
    Abstract base class for marker item headers.
    """

    pass


class MarkerItemHeaderPre4(MarkerItemHeader):
    """
    Marker Items for files in Version 3, very likely down to version 2.
    """

    is_valid_for_revision = pre_4()

    _versioned_fields = [
        ("lSample", ctypes.c_int32, V_20a),
        ("fSelected", ctypes.c_int16, V_35x),
        ("fTextLocked", ctypes.c_int16, V_20a),
        ("fPositionLocked", ctypes.c_int16, V_20a),
        ("nTextLength", ctypes.c_int16, V_20a),
    ]

    @property
    def text_length(self):
        if self.file_revision < V_35x:
            return self._struct.nTextLength
        else:
            return self._struct.nTextLength + 1

    @property
    def effective_len_bytes(self):
        return self.struct_length + self.text_length

    @property
    def sample_index(self):
        return self._struct.lSample

    @property
    def channel_number(self):
        """None means it's a global marker"""
        return None

    @property
    def date_created_ms(self):
        """markers don't get "Date created" until v. 4.4.0"""
        return None

    @property
    def type_code(self):
        """These markers don't get type_codes, but it's OK"""
        return None


class MarkerItemHeaderPost4(MarkerItemHeader):
    """
    Marker Items for files in Version 4 onwards.
    """

    is_valid_for_revision = post_4()

    # Define the structure fields
    _versioned_fields = [
        ("lSample", ctypes.c_uint32, V_400B),
        ("Unknown", ctypes.c_byte * 4, V_400B),
        ("nChannel", ctypes.c_int16, V_400B),
        ("sMarkerStyle", ctypes.c_char * 4, V_400B),
        ("llDateCreated", ctypes.c_uint64, V_440),
        ("Unknown3", ctypes.c_byte * 8, V_42x),
        ("nTextLength", ctypes.c_int16, V_400B),
    ]

    @property
    def text_length(self):
        """Get the text length, including the trailing null."""
        return self._struct.nTextLength

    @property
    def effective_len_bytes(self):
        return self.struct_length + self.text_length

    @property
    def sample_index(self):
        """Get the sample index."""
        return self._struct.lSample

    @property
    def channel_number(self):
        """Get the channel number, or None if it's a global marker."""
        chan = self._struct.nChannel
        return None if chan == -1 else chan

    @property
    def date_created_ms(self):
        """Get the date when the marker was created (in ms since 1970-01-01)."""
        if self.file_revision < V_440:
            return None
        return self._struct.llDateCreated

    @property
    def type_code(self):
        """Get the type code, decoded as a string."""
        return self._struct.sMarkerStyle.decode(self.encoding, errors="ignore")
