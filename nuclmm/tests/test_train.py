# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of nuclmm (http://github.com/standage/nuclmm) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function, division
import nuclmm
from nuclmm import (parse_fasta, data_file)
from nuclmm.markovchain import (MarkovChain, MarkovChainOrderMismatch,
                                MarkovChainNormalizationError)
import pytest


def test_order_5():
    model = MarkovChain(order=5)
    with nuclmm.open(data_file('bogus.fa.gz'), 'r') as inseq:
        for defline, sequence in nuclmm.parse_fasta(inseq):
            model.train(sequence)
    model.normalize()

    testmodel = MarkovChain(order=5)
    with nuclmm.open(data_file('bogus.order5.mm'), 'r') as infile:
        testmodel.load(infile)
    assert model == testmodel
