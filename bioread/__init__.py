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

from bioread import reader

from ._metadata import version as __version__, author as __author__  # noqa


def read(filelike, channel_indexes=None):
    """
    Read a file (either an IO object or a filename) and return a Datafile.

    channel_indexes:    A list of integer channel numbers. Other channels will
                        have empty data.
    target_chunk_size:  A guide for the number of bytes to read at a time.
    """
    return reader.Reader.read(filelike, channel_indexes).datafile

# Deprecated; provided for compatibility with previous versions.
read_file = read


def read_headers(filelike):
    """
    Read only the headers of a file, returns a Datafile with empty channels.
    """
    return reader.Reader.read_headers(filelike).datafile


def reader_for_streaming(io):
    """
    Read the headers of a file, return a Reader object that will allow you to
    stream the data in chunks with stream().
    """
    if not hasattr(io, 'read'):
        raise TypeError('{0} must be an opened file.'.format(io))
    if hasattr(io, 'encoding'):
        raise TypeError('{0} must be opened in binary mode'.format(io))
    return reader.Reader.read_headers(io)
