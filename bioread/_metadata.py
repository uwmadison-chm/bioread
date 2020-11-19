# coding: utf8
# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2020 Board of Regents of the University of Wisconsin System
#
# Written Nate Vack <njvack@wisc.edu> with research from John Ollinger
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison
# Project home: http://github.com/njvack/bioread

# NOTE: This file must not import anything, or it will break installation.

version_tuple = (2, 1, 0)
version = ".".join([str(p) for p in version_tuple])
version_description = "bioread {0}".format(version)

author = "Nate Vack"
author_email = "njvack@wisc.edu"
license = "GPL 2.0"
copyright = "Copyright 2016 Boards of Regent of the University of Wisconsin System"  # noqa
url = "https://github.com/uwmadison-chm/bioread"
