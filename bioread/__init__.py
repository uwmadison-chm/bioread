# coding: utf8
# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2025 Board of Regents of the University of Wisconsin System
#
# Written Nate Vack <njvack@wisc.edu> with research from John Ollinger
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison

from typing import Union, List, Optional, BinaryIO
import os

from bioread import reader
from bioread.biopac import Datafile

from ._metadata import __version__, author as __author__  # noqa


def read(
    filelike: Union[str, os.PathLike, BinaryIO],
    channel_indexes: Optional[List[int]] = None,
) -> Datafile:
    """
    Read a file (either an IO object or a filename) and return a Datafile.

    Parameters
    ----------
    filelike : str, PathLike, or file-like object
        The file to read
    channel_indexes : list of int, optional
        A list of integer channel numbers. Other channels will have empty data.

    Returns
    -------
    Datafile
        The read datafile
    """
    return reader.Reader.read(filelike, channel_indexes).datafile


# Deprecated; provided for compatibility with previous versions.
read_file = read


def read_headers(filelike: Union[str, os.PathLike, BinaryIO]) -> Datafile:
    """
    Read only the headers of a file, returns a Datafile with empty channels.

    Parameters
    ----------
    filelike : str, PathLike, or file-like object
        The file to read

    Returns
    -------
    Datafile
        The datafile with headers read
    """
    return reader.Reader.read_headers(filelike).datafile


def reader_for_streaming(io: BinaryIO) -> reader.Reader:
    """
    Read the headers of a file, return a Reader object that will allow you to
    stream the data in chunks with stream().

    Parameters
    ----------
    io : file-like object
        An opened file in binary mode

    Returns
    -------
    Reader
        A reader object for streaming data

    Raises
    ------
    TypeError
        If io is not a file-like object or is not opened in binary mode
    """
    if not hasattr(io, "read"):
        raise TypeError(f"{io} must be an opened file.")
    if hasattr(io, "encoding"):
        raise TypeError(f"{io} must be opened in binary mode")
    return reader.Reader.read_headers(io)
