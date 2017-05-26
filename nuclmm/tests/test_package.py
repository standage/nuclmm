# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of nuclmm (http://github.com/standage/nuclmm) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function
import pytest
import nuclmm


def test_open_badmode():
    with pytest.raises(ValueError) as ve:
        nuclmm.open('bogusfile', 'z')
    assert 'invalid mode "z"' in str(ve)


def test_fasta():
    infilename = nuclmm.data_file('threeseqs.fa.gz')
    infile = nuclmm.open(infilename, 'r')
    seqs = [seq for defline, seq in nuclmm.parse_fasta(infile)]
    assert len(seqs) == 3
