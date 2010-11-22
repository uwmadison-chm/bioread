# coding: utf8
# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2010 Board of Regents of the University of Wisconsin System
#
# Written by John Ollinger <ollinger@wisc.edu> and Nate Vack <njvack@wisc.edu>
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison
# Project home: http://github.com/njvack/bioread

# Some basic useful stuff that doesn't fit elegantly elsewhere.

# Simple LCM routine adapted from:
# http://stackoverflow.com/questions/147515/least-common-multiple-for-3-or-more-numbers

def gcd(a, b):
    while b:
        a, b = b, a % b
    return a

def lcm(*ar):
    if len(ar) > 2:
        return lcm(ar[0], lcm(*ar[1:]))
    elif len(ar) == 2:
        return (ar[0] * ar[1]) // gcd(ar[0], ar[1])
    else:
        return ar[0] 