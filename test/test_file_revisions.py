# coding: utf8
# Part of the bioread package for reading BIOPAC data.
#
# Copyright (c) 2025 Board of Regents of the University of Wisconsin System
#
# Written Nate Vack <njvack@wisc.edu> with research from John Ollinger
# at the Waisman Laboratory for Brain Imaging and Behavior, University of
# Wisconsin-Madison

# Tests for the file_revisions functions


from bioread import file_revisions as rev


def test_version_stringify_examples():
    assert "5.0.1" == rev.stringify_version("V_501")
    assert "unknown early version" == rev.stringify_version("V_ALL")
    assert "4.0.0.B" == rev.stringify_version("V_400B")


def test_version_string_guesses():
    assert "5.0.1" == rev.version_string_guess(132)
    assert "after 5.0.1" == rev.version_string_guess(500)
    assert "unknown early version" == rev.version_string_guess(-1)
    assert "between 4.4.0 and 5.0.1" == rev.version_string_guess(131)
    assert "between unknown early version and 2.0.a" == rev.version_string_guess(1)
