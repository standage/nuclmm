#!/usr/bin/env python
from __future__ import print_function, division
from nuclmm import (parse_fasta, MarkovChain)
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


def test_init_probs():
    model = MarkovChain(order=2)
    with open('aaat.fa', 'r') as infile:
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
