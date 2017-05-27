# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of nuclmm (http://github.com/standage/nuclmm) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function, division
from io import StringIO
import pytest
import nuclmm
from nuclmm import (parse_fasta, data_file)
from nuclmm.markovchain import (MarkovChain, MarkovChainOrderMismatch,
                                MarkovChainNormalizationError, draw, nuclpairs)


def test_basic():
    assert MarkovChain().order == 1
    assert MarkovChain(order=5).order == 5

    with pytest.raises(AssertionError):
        model = MarkovChain(order=0)

    with pytest.raises(AssertionError):
        model = MarkovChain(order=-1)

    model = MarkovChain(order=3)
    assert str(model) == '[\n    {},\n    {}\n]'


def test_draw():
    assert draw({}) is None


def test_nucl_pairs():
    assert list(nuclpairs('A')) == list()
    assert list(nuclpairs('ACGTTTT')) == [('A', 'C'), ('C', 'G'), ('G', 'T'),
                                          ('T', 'T'), ('T', 'T'), ('T', 'T')]


def test_train():
    model = MarkovChain(order=2)
    model.train('AAATAAATAAAT')
    assert sorted(model._transitions) == ['AA', 'AT', 'TA']
    assert list(model._initial_states) == ['AA']

    model = MarkovChain(order=6)
    model.train('TNTAAAGGGGCTTTCGAGTGTGATCAAACGGCTCCAAACCCACTAAACAAAATACGCTCTG'
                'TATTGGGCCCAAAGATTTCGCGTGTTGATCGACGGGCGATCCGTTCGCACGCA')

    model = MarkovChain(order=6)
    model.train('TCTAAAGGGGCTTTCGAGTGTGATCAAACGGCTCCAAACCCACTAAACAAAATACGCTCTG'
                'TATTGGGCCCAAAGATTTCGCGTNTTGATCGACGGGCGATCCGTTCGCACGCA')


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


def test_inequality():
    infile1 = nuclmm.open(data_file('bogus.order2.mm'), 'r')
    infile2 = nuclmm.open(data_file('bogus.order2.difftrans.mm'), 'r')
    infile3 = nuclmm.open(data_file('bogus.order2.diffprob.mm'), 'r')
    infile4 = nuclmm.open(data_file('bogus.order5.mm'), 'r')

    model1 = MarkovChain(order=2)
    model2 = MarkovChain(order=2)
    model3 = MarkovChain(order=2)
    model4 = MarkovChain(order=5)

    model1.load(infile1)
    model2.load(infile2)
    model3.load(infile3)
    model4.load(infile4)

    assert model1 != model2
    assert model1 != model3
    assert model1 != model4


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


def test_open():
    model = MarkovChain(order=3)
    infile = data_file('aaat.mm')
    model.open(infile)

    testmodel = MarkovChain(order=3)
    with nuclmm.open(infile, 'r') as infile:
        testmodel.load(infile)

    assert model == testmodel


def test_save():
    model = MarkovChain(order=3)
    model.train('AAATAAATAAAT')
    iostream = StringIO()
    model.save(iostream)

    iostream.seek(0)
    loadmodel = MarkovChain(order=3)
    loadmodel.load(iostream)

    testmodel = MarkovChain(order=3)
    testmodel.load(nuclmm.open(data_file('aaat.mm'), 'r'))

    assert loadmodel == testmodel


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


def test_simulate():
    model = MarkovChain(order=5)
    with nuclmm.open(data_file('bogus.order5.mm'), 'r') as infile:
        model.load(infile)

    seqs = [s for s in model.simulate(100, 100, ignore_inits=True)]
    assert len(seqs) == 100
    for seq in seqs:
        assert len(seq) == 100
