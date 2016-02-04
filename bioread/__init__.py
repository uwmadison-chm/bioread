# coding: utf8
# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2016 Board of Regents of the University of Wisconsin System
#
# Written Nate Vack <njvack@wisc.edu> with research from John Ollinger
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison
# Project home: http://github.com/njvack/bioread

from __future__ import absolute_import

from bioread.readers import AcqReader

from ._metadata import version as __version__, author as __author__  # noqa


def read_file(filelike):
    """
    Read a file (either an IO object or a filename) and return a
    biopac.Datafile object. Simply a shorthand for
    bioread.readers.AcqReader.read_file()
    """
    return AcqReader.read_file(filelike)
