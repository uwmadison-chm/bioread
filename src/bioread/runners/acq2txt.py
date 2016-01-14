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

# This contains the entry point for an executable to convert BIOPAC
# AcqKnowledge files into text files
""" Write the data from an AcqKnowledge file channel to a text file.

Usage:
  acq2txt [options] <acq_file>
  acq2txt -h | --help
  acq2txt --version

Options:
  --version          show program's version number and exit
  -h, --help         show this help message and exit
  --channel=CHANNEL  channel number to extract [default: 0]

Writing more than one channel is not supported at the current time, because
different channels can have different sampling rates, and it's hard to know
what to do in those cases.
"""

import sys

from bioread.vendor.docopt import docopt

import bioread
from bioread.writers import TxtWriter
from bioread import version


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

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
        pargs = docopt(
            __doc__,
            self.argv,
            version=version.description)
        infile = pargs['<acq_file>']
        try:
            data = bioread.read_file(infile)
        except:
            sys.stderr.write("Error reading %s\n" % infile)
            sys.exit(1)
        try:
            chan = data.channels[int(pargs['--channel'])]
        except:
            sys.stderr.write(
                "Channel %s out of bounds -- max: %s\n" %
                (pargs['--channel'], len(data.channels) - 1))
            sys.exit(2)
        try:
            TxtWriter.write_file(chan, sys.stdout)
        except:
            sys.stderr.write("Notice: Not all data written\n")
            sys.exit(-1)

        sys.stderr = old_err


if __name__ == '__main__':
    main()
