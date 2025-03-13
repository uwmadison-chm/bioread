# coding: utf8
# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2025 Board of Regents of the University of Wisconsin System
#
# Written Nate Vack <njvack@wisc.edu> with research from John Ollinger
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison

# This is an enumeration of known file revisions, for use in header parsing.

V_ALL = 0
V_20a = 30
V_20b = 31
V_20r = 32
V_207 = 33
V_30r = 34
V_303 = 35
V_35x = 36
V_36x = 37
V_370 = 38
V_373 = 39
V_381 = 41
V_37P = 42
V_382 = 43
V_38P = 44
V_390 = 45
V_400B = 61  # Earliest 'late version' one I've seen so far...
V_400 = 68
V_401 = 76
V_402 = 78  # Unsure about exact versions... but this is Post-Some-Big-Change.
V_41a = 80  # But I'm honestly guessing about some of these.
V_410 = 83
V_411 = 84
V_420 = 108
V_42x = 121  # Not yet sure what version this is
V_430 = 124
V_440 = 128
V_501 = 132

# Make a dictionary mapping revison numbers to the version string.
# This is a silly amount of code for a very marginal need


def stringify_version(version_name):
    """
    Convert something like V_400 to "4.0.0"
    """
    no_v = version_name.replace("V_", "")
    if no_v.startswith("A"):
        return "unknown early version"
    return ".".join(no_v)


# wooooo looping over globals(), I kind of hate this

REVISION_TO_VERSION = {
    revision: stringify_version(ver_const_name)
    for ver_const_name, revision in globals().items()
    if isinstance(revision, int) and ver_const_name.startswith("V_")
}


def version_string_guess(revision):
    """
    If we have revision number in our dictionary, return the version string.
    Otherwise, if the revision number is less than the first version, return
    "before {version}". If it's between two versions, return "between {version1} and {version2}".
    Otherwise, return "after {version}".
    """
    sorted_revisions = sorted(REVISION_TO_VERSION.keys())  # int list
    if revision in REVISION_TO_VERSION:
        return REVISION_TO_VERSION[revision]
    if revision < sorted_revisions[0]:
        return f"unknown early version"
    if revision > sorted_revisions[-1]:
        return f"after {REVISION_TO_VERSION[sorted_revisions[-1]]}"
    for i, check_revision in enumerate(sorted_revisions):
        if revision < check_revision:
            print(i)
            low_bound = REVISION_TO_VERSION[sorted_revisions[i - 1]]
            high_bound = REVISION_TO_VERSION[sorted_revisions[i]]
            return f"between {low_bound} and {high_bound}"
