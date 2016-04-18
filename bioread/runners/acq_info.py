#!/usr/bin/env python
# coding: utf8
# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2016 Board of Regents of the University of Wisconsin System
#
# Written Nate Vack <njvack@wisc.edu> with research from John Ollinger
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison
# Project home: http://github.com/njvack/bioread

# This contains the entry point for an executable to print basic information
# about an AcqKnowledge file.

"""Print some information about an AcqKnowledge file.

Usage:
    acq_info [options] <acq_file>
    acq_info -h | --help
    acq_info --version

Options:
  -d, --debug  print lots of debugging data

Note: Using - for <acq_file> reads from stdin.

"""

import sys

from bioread.vendor.six import BytesIO
from bioread.vendor.docopt import docopt

from bioread.reader import Reader
from bioread import version


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    air = AcqInfoRunner(argv)
    air.run()


class AcqInfoRunner(object):

    def __init__(self, argv, out=None, err=None):
        self.argv = argv
        if out is None:
            out = sys.stdout
        self.out = out
        if err is None:
            err = sys.stderr
        self.err = err

    def run(self):
        old_out = sys.stdout
        old_err = sys.stderr
        sys.stdout = self.out
        sys.stderr = self.err

        pargs = docopt(
            __doc__,
            self.argv,
            version=version.description)

        df = None
        infile = pargs['<acq_file>']
        try:
            if infile == '-':
                df = BytesIO(sys.stdin.read())
            else:
                df = open(infile, 'rb')
        except:
            sys.stderr.write("Error reading {0}\n".format(infile))
            sys.exit(1)

        self.reader = Reader(df)
        try:
            self.reader._read_headers()
        except:
            sys.stderr.write("Error reading headers!\n")
            # Don't exit here; it'll still print what it can.

        if pargs['--debug']:
            self.__print_debug()
        else:
            self.__print_simple()

        sys.stderr = old_err
        sys.stdout = old_out

    def __print_simple(self):
        gh = self.reader.graph_header
        chs = self.reader.channel_headers
        cdhs = self.reader.channel_dtype_headers
        print("File revision: %s" % gh.file_revision)
        print("Sample time: %s" % gh.sample_time)
        print("Compressed: %s" % gh.compressed)
        print("Number of channels: %s" % gh.channel_count)
        for ch, cdh in zip(chs, cdhs):
            print("%s:" % ch.name)
            print("\tUnits: %s" % ch.units)
            print("\tNumber of samples: %s" % ch.point_count)
            print("\tFrequency divider: %s" % ch.frequency_divider)

    def __print_debug(self):
        gh = self.reader.graph_header
        fh = self.reader.foreign_header
        chs = self.reader.channel_headers
        cdhs = self.reader.channel_dtype_headers

        print("Graph header starts at offset %s" % gh.offset)
        print(gh.data)
        print("\n")
        for i, ch in enumerate(chs):
            print("Channel header %s starts at offset %s" % (i, ch.offset))
            print(ch.data)
            print("\n")
        print("Foreign header starts at offset %s" % fh.offset)
        print(fh.data)
        print("\n")
        for i, cdh in enumerate(cdhs):
            print(
                "Channel dtype header %s starts at offset %s" %
                (i, cdh.offset))
            print(cdh.data)
            print("\n")

        if not gh.compressed:
            print("Data starts at offset %s" % self.reader.data_start_offset)
        else:
            mch = self.reader.main_compression_header
            cchs = self.reader.channel_compression_headers
            print("Main compression header starts at offset %s" % mch.offset)
            print("\n")
            for i, cch in enumerate(cchs):
                print(
                    "Channel compression header %s starts at offset %s" %
                    (i, cch.offset))
                print(cch.data)
                print(
                    "Compressed data starts at offset %s" %
                    cch.compressed_data_offset)
                print("\n")


if __name__ == '__main__':
    main()
