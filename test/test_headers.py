# coding: utf8
# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2025 Board of Regents of the University of Wisconsin System
#
# Written by Nate Vack <njvack@wisc.edu> with research from John Ollinger
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison, and later the Center for Healthy Minds

import bioread.file_revisions as rev
from bioread import headers as bh


def test_struct_class_is_shared_between_identical_headers():
    """Building the ctypes class is expensive; we should only do it once."""
    first = bh.ChannelHeader.for_revision(rev.V_400B, "<")
    second = bh.ChannelHeader.for_revision(rev.V_400B, "<")
    assert first is not second
    assert first._struct_class is second._struct_class


def test_struct_class_differs_by_byte_order():
    little = bh.ChannelHeader.for_revision(rev.V_400B, "<")
    big = bh.ChannelHeader.for_revision(rev.V_400B, ">")
    assert little._struct_class is not big._struct_class


def test_struct_class_differs_by_revision():
    """V_373 adds vertPrecision to the pre-4 channel header."""
    older = bh.ChannelHeader.for_revision(rev.V_370, "<")
    newer = bh.ChannelHeader.for_revision(rev.V_373, "<")
    assert older._struct_class is not newer._struct_class

    older_fields = [f[0] for f in older._struct_class._fields_]
    newer_fields = [f[0] for f in newer._struct_class._fields_]
    assert "vertPrecision" not in older_fields
    assert "vertPrecision" in newer_fields
    assert newer.struct_length > older.struct_length


def test_struct_class_differs_by_header_class():
    channel = bh.ChannelHeader.for_revision(rev.V_400B, "<")
    marker = bh.MarkerHeader.for_revision(rev.V_400B, "<")
    assert channel._struct_class is not marker._struct_class


def test_cached_struct_class_still_unpacks_correctly():
    """A shared struct class must not mean shared parsed data."""
    header_class = bh.MarkerHeaderPre4
    length = bh.MarkerHeader.for_revision(rev.V_20a, "<").struct_length

    first = header_class(rev.V_20a, "<", "utf-8", _from_factory=True)
    second = header_class(rev.V_20a, "<", "utf-8", _from_factory=True)

    first.unpack_from_bytes(b"\x01\x00\x00\x00\x02\x00\x00\x00"[:length])
    second.unpack_from_bytes(b"\x03\x00\x00\x00\x04\x00\x00\x00"[:length])

    assert first._struct_class is second._struct_class
    assert first.marker_count == 2
    assert second.marker_count == 4


def test_byte_order_is_respected_by_cached_class():
    length = bh.MarkerHeader.for_revision(rev.V_20a, "<").struct_length
    raw = b"\x00\x00\x00\x01\x00\x00\x00\x02"[:length]

    little = bh.MarkerHeaderPre4(rev.V_20a, "<", "utf-8", _from_factory=True)
    big = bh.MarkerHeaderPre4(rev.V_20a, ">", "utf-8", _from_factory=True)
    little.unpack_from_bytes(raw)
    big.unpack_from_bytes(raw)

    assert big.marker_count == 2
    assert little.marker_count == 0x02000000
