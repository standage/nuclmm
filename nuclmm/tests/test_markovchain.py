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


def test_basic():
    assert MarkovChain().order == 1
    assert MarkovChain(order=5).order == 5

    with pytest.raises(AssertionError):
        model = MarkovChain(order=0)

    with pytest.raises(AssertionError):
        model = MarkovChain(order=-1)

    model = MarkovChain(order=3)
    assert str(model) == '[\n    {},\n    {}\n]'


def test_train():
    model = MarkovChain(order=2)
    model.train('AAATAAATAAAT')
    assert sorted(model._transitions) == ['AA', 'AT', 'TA']
    assert list(model._initial_states) == ['AA']


def test_equality():
    assert MarkovChain(order=4) != MarkovChain(order=2)

    model1 = MarkovChain(order=4)
    model1.train('AAGCGTCCGGTGTATTGTCTAGCTGGCTAAATTTTCATACTAA')
    model1.train('TCTACCTAACGCCTTAACGATATGCGAAGCGACTACTCTCACT')

    model2 = MarkovChain(order=4)
    model2.train('AAGCGTCCGGTGTATTGTCTAGCTGGCTAAATTTTCATACTAA')
    model2.train('TCTACCTAACGCCTTAACGATATGCGAAGCGACTACTCTCACT')

    model3 = MarkovChain(order=4)
    model3.train('CATACCCCGTGTACAGCGTTATAACCATTCGCTAAAAGTGACA')
    model3.train('GGGTGAACTTACGTGATATGAGCCTACTGCAGGGTCGCCTACG')

    assert model1 == model2
    assert model1 != model3


def test_init_probs():
    model = MarkovChain(order=2)
    infile = data_file('aaat.fa')
    with nuclmm.open(infile, 'r') as infile:
        for defline, sequence in parse_fasta(infile):
            model.train(sequence)
    state, prob = next(model.initial_probs)
    assert state == 'AA'
    assert prob == 1.0


def test_trans_probs():
    model = MarkovChain(order=2)
    model.train('AAATAAATAAAT')
    for state, nextnucl, prob in model.transition_probs:
        if state == 'AA' and nextnucl in ['A', 'T']:
            assert prob == 0.5, (state, nextnucl, prob)


def test_load():
    model = MarkovChain(order=3)
    infile = data_file('aaat.mm')
    with nuclmm.open(infile, 'r') as infile:
        model.load(infile)
    assert sorted(model._initial_states.keys()) == ['AAA']
    assert sorted(model._transitions.keys()) == ['AAA', 'AAT', 'ATA', 'TAA']


def test_order_mismatch():
    model = MarkovChain(order=2)
    infile = data_file('aaat.mm')
    with pytest.raises(MarkovChainOrderMismatch) as mcom:
        with nuclmm.open(infile, 'r') as infile:
            model.load(infile)
    assert 'state has length 3 but model declared 2-order' in str(mcom)


def test_init_state_norm():
    model = MarkovChain(order=4)
    infile = data_file('init-state-norm.mm')
    with pytest.raises(MarkovChainNormalizationError) as mvne:
        with nuclmm.open(infile, 'r') as infile:
            model.load(infile)
    assert 'initial state probabilites should sum to 1.0' in str(mvne)


def test_state_trans_norm():
    model = MarkovChain(order=4)
    infile = data_file('state-trans-norm.mm')
    with pytest.raises(MarkovChainNormalizationError) as mvne:
        with nuclmm.open(infile, 'r') as infile:
            model.load(infile)
    assert 'transition probabilities for TGAA should sum to 1.0' in str(mvne)
