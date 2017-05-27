# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of nuclmm (http://github.com/standage/nuclmm) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function
import pytest
import subprocess
import nuclmm
from nuclmm import data_file


def test_open_badmode():
    with pytest.raises(ValueError) as ve:
        nuclmm.open('bogusfile', 'z')
    assert 'invalid mode "z"' in str(ve)


def test_fasta():
    infilename = nuclmm.data_file('threeseqs.fa.gz')
    infile = nuclmm.open(infilename, 'r')
    seqs = [seq for defline, seq in nuclmm.parse_fasta(infile)]
    assert len(seqs) == 3


def test_main(capfd):
    traincmd = './cli train --order 4 {}'.format(data_file('aaat.fa'))
    trainargs = traincmd.split()
    trainproc = subprocess.Popen(trainargs, stdout=subprocess.PIPE)

    simcmd = './cli simulate --order 4 --numseqs 1 --seqlen 32 --seed 42 -'
    simargs = simcmd.split()
    simproc = subprocess.Popen(simargs, stdin=trainproc.stdout)

    simproc.communicate()
    out, err = capfd.readouterr()
    assert 'AAATAAATAAATAAATAAATAAATAAATAAAT' in out
