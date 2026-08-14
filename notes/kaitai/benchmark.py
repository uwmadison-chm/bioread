#!/usr/bin/env python
# coding: utf8
# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2025 Board of Regents of the University of Wisconsin System

"""Compare Kaitai and ctypes header parsing on cost and on bytes touched.

Usage:
    python notes/kaitai/benchmark.py
"""

import ctypes
import io
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO_ROOT = HERE.parents[1]
sys.path.insert(0, str(REPO_ROOT))
sys.path.insert(0, str(HERE))

from kaitaistruct import KaitaiStream  # noqa: E402

from bioread import headers as bh  # noqa: E402
from bioread.header_reader import GraphHeaderReader  # noqa: E402

try:
    from acq import Acq
except ImportError:
    sys.exit("No generated parser found -- see validate.py for the build steps.")

SAMPLE = REPO_ROOT / "test" / "data" / "physio" / "physio-4.4.0.acq"
REVISION = 128  # what physio-4.4.0.acq reports
REPEATS = 500


class CountingFile(io.RawIOBase):
    """A read-only file wrapper that tallies how much actually gets read."""

    def __init__(self, path):
        self.f = open(path, "rb")
        self.bytes_read = 0
        self.reads = 0
        self.seeks = 0

    def read(self, n=-1):
        data = self.f.read(n)
        self.bytes_read += len(data)
        self.reads += 1
        return data

    def readinto(self, b):
        data = self.f.read(len(b))
        b[: len(data)] = data
        self.bytes_read += len(data)
        self.reads += 1
        return len(data)

    def seek(self, offset, whence=0):
        self.seeks += 1
        return self.f.seek(offset, whence)

    def tell(self):
        return self.f.tell()

    def seekable(self):
        return True

    def readable(self):
        return True


def timeit(label, fn, repeats=REPEATS):
    fn()  # warm up
    start = time.perf_counter()
    for _ in range(repeats):
        fn()
    per_call = (time.perf_counter() - start) / repeats
    print(f"  {label:<44} {per_call * 1e6:8.1f} us")
    return per_call


def report_io():
    print("Bytes touched to read the graph header:\n")

    counting = CountingFile(SAMPLE)
    parsed = Acq(KaitaiStream(counting))
    _ = parsed.graph_header.n_channels
    print(
        f"  kaitai   {counting.bytes_read:>7,} bytes  "
        f"{counting.reads:>5} reads  {counting.seeks:>3} seeks"
    )

    counting = CountingFile(SAMPLE)
    header = GraphHeaderReader.bootstrap(counting)
    _ = header.channel_count
    print(
        f"  ctypes   {counting.bytes_read:>7,} bytes  "
        f"{counting.reads:>5} reads  {counting.seeks:>3} seeks"
    )


def report_speed():
    raw = SAMPLE.read_bytes()

    print("\nParsing one graph header:\n")

    def ctypes_everything():
        header = bh.GraphHeader.for_revision(REVISION, "<", "utf-8")
        header.unpack_from_bytes(raw[: header.struct_length])
        return header.channel_count

    prebuilt = bh.GraphHeader.for_revision(REVISION, "<", "utf-8")._struct_class
    buffer = raw[: ctypes.sizeof(prebuilt)]

    def ctypes_build_class():
        return bh.GraphHeader.for_revision(REVISION, "<", "utf-8")._struct_class

    def ctypes_parse_only():
        return prebuilt.from_buffer_copy(buffer).nChannels

    def kaitai_parse():
        return Acq(KaitaiStream(io.BytesIO(raw))).graph_header.n_channels

    timeit("ctypes: for_revision + class + parse", ctypes_everything)
    timeit("ctypes: build the struct class only", ctypes_build_class)
    timeit("ctypes: parse only, class reused", ctypes_parse_only)
    timeit("kaitai: construct + parse", kaitai_parse)


def report_marker_scaling(count=10000):
    print(f"\nReading {count:,} marker item headers:\n")

    blank = bytes(64)

    def via_reader():
        for _ in range(count):
            header = bh.MarkerItemHeader.for_revision(REVISION, "<", "utf-8")
            header.unpack_from_bytes(blank[: header.struct_length])

    prebuilt = bh.MarkerItemHeader.for_revision(REVISION, "<", "utf-8")._struct_class
    buffer = blank[: ctypes.sizeof(prebuilt)]

    def class_reused():
        for _ in range(count):
            prebuilt.from_buffer_copy(buffer)

    timeit("bioread MarkerItemHeader.for_revision", via_reader, repeats=3)
    timeit("ctypes with the class reused", class_reused, repeats=3)


if __name__ == "__main__":
    print(f"sample: {SAMPLE.name} ({SAMPLE.stat().st_size:,} bytes)\n")
    report_io()
    report_speed()
    report_marker_scaling()
