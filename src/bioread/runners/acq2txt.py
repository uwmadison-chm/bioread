#!/usr/bin/env python
# coding: utf8
# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2010 Board of Regents of the University of Wisconsin System
#
# Written by John Ollinger <ollinger@wisc.edu> and Nate Vack <njvack@wisc.edu>
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison
# Project home: http://github.com/njvack/bioread

# This contains the entry point for an executable to convert BIOPAC
# AcqKnowledge files into Matlab files.

import sys
import os.path
import StringIO
from optparse import OptionParser

import bioread
from bioread.writers import TxtWriter
from bioread.version import version_str


def main(argv = None):
    if argv is None:
        argv = sys.argv

    amr = AcqToTxtRunner(argv)
    amr.run()


class AcqToTxtRunner(object):
    """The little wrapper class that converts acq files to mat files"""

    def __init__(self, argv, err=None):
        self.argv = argv
        if err is None:
            err = sys.stderr
        self.err = err

    def run(self):
        old_err = sys.stderr
        sys.stderr = self.err
        self.parser = self.__make_parser()
        options, args = self.parser.parse_args(self.argv[1:])
        if len(args) <> 1:
            self.parser.error(
                "Must specify ACQ_FILE.\n"+
                "Try --help for more instructions.")
        try:
            infile = args[0]
            if infile == '-':
                infile = StringIO.StringIO(sys.stdin.read())
            data = bioread.read_file(infile)
        except:
            sys.stderr.write("Error reading %s\n" % args[0])
            sys.exit(1)
        try:
            chan = data.channels[options.channel]
        except:
            sys.stderr.write("Channel %s out of bounds -- max: %s\n" %
                (options.channel, len(data.channels)))
            sys.exit(2)
        try:
            TxtWriter.write_file(chan, sys.stdout)
        except:
            sys.stderr.write("Notice: Not all data written\n")
            sys.exit(-1)

        sys.stderr = old_err

    def __make_parser(self):
        parser = OptionParser(
            "Usage: %prog [options] ACQ_FILE",
            version="bioread %s" % version_str(),
            epilog="Note: Using - for ACQ_FILE reads from stdin.")
        parser.add_option("--channel", dest="channel", default=0,
            action="store", help="channel number to extract (default=0)",
            type="int")

        return parser

if __name__ == '__main__':
    main()
