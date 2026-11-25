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
    assert rev.stringify_version("V_501") == "5.0.1"
    assert rev.stringify_version("V_ALL") == "unknown early version"
    assert rev.stringify_version("V_400B") == "4.0.0.B"


def test_version_string_guesses():
    assert rev.version_string_guess(rev.V_501) == "AcqKnowledge 5.0.1"
    assert rev.version_string_guess(rev.V_36x) == "BSL PRO 3.6.x"
    assert (
        rev.version_string_guess(rev.V_501 + 1) == "Probably AcqKnowledge, after 5.0.1"
    )
    assert (
        rev.version_string_guess(rev.V_ALL - 1)
        == "Probably AcqKnowledge, unknown early version"
    )
    assert (
        rev.version_string_guess(rev.V_440 + 1)
        == "Probably AcqKnowledge, between 4.4.0 and 5.0.1"
    )
    assert (
        rev.version_string_guess(rev.V_ALL + 1)
        == "Probably AcqKnowledge, between unknown early version and 2.0.a"
    )
