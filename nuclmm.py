#!/usr/bin/env python
from __future__ import print_function, division
from collections import deque, defaultdict
from itertools import chain
import argparse
import re
import sys
import pytest


# -----------------------------------------------------------------------------
# Utility functions
# -----------------------------------------------------------------------------

def parse_fasta(data):
    """Stolen shamelessly from http://stackoverflow.com/a/7655072/459780."""
    name, seq = None, []
    for line in data:
        line = line.rstrip()
        if line.startswith('>'):
            if name:
                yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name:
        yield (name, ''.join(seq))


def nuclpairs(sequence):
    if len(sequence) < 2:
        return None
    prev = sequence[0]
    for nucl in sequence[1:]:
        yield prev, nucl
        prev = nucl


def nucldict():
    return {n: 0 for n in 'ACGT'}


# -----------------------------------------------------------------------------
# Class definition
# -----------------------------------------------------------------------------

class NuclCompModel(object):
    """
    Class for training Markov models of nucleotide composition.

    In this context an Nth order Markov model is defined as follows: first, a
    set of *transition probabilities*, i.e., the probability that the next
    nucleotide will be *X* given that the current state (*N* nucleotides) is
    *Y*; and second, a set of *initial state probabilities*, i.e., the
    probability of observing *Z* as the first *N* nucleotides of a given
    sequence.

    This class supports models of arbitrary order (N >= 1).

    Probabilities are computed by observing invoking the **consume** command on
    one or more training sequences.
    """

    def __init__(self, order=1):
        assert order > 0
        self._initial_states = defaultdict(int)
        self._transitions = defaultdict(nucldict)
        self._nuclstate = deque()
        self._order = int(order)

    def train(self, sequence):
        """
        Train the model on the given sequence.

        :param sequence: string-like object containing a nucleotide sequence
        """
        self._record_initial_state(sequence)
        for subseq in re.split('[^ACGT]+', sequence.upper()):
            for nucl, nextnucl in nuclpairs(subseq):
                self._nuclstate.append(nucl)
                if len(self._nuclstate) < self._order:
                    continue
                state = ''.join([n for n in self._nuclstate])
                self._transitions[state][nextnucl] += 1
                self._nuclstate.popleft()
            self._nuclstate = deque()  # Reset for next sequence/subsequence

    def load(self, instream, num_header_lines=1):
        for _ in range(num_header_lines):
            _ = next(instream)
        for line in instream:
            if line.strip() == '':
                continue
            state, nextnucl, prob = line.strip().split(',')
            assert len(state) == self.order
            prob = float(prob)
            if nextnucl == '':
                self._initial_states[state] = prob
            else:
                self._transitions[state][nextnucl] = prob

    def _record_initial_state(self, sequence):
        init = sequence[:self._order]
        if re.search('[^ACGT]', init):
            message = ('Initial state "{}" includes an ambiguous/erroneous'
                       'nucleotide, ignoring'.format(init))
            print(message, file=sys.stderr)
        else:
            self._initial_states[init] += 1

    @property
    def order(self):
        return self._order

    @property
    def initial_probs(self):
        norm = sum(self._initial_states.values())
        for initseq in sorted(self._initial_states):
            count = self._initial_states[initseq]
            prob = count / norm
            yield initseq, prob

    @property
    def transition_probs(self):
        for stateseq in sorted(self._transitions):
            counts = self._transitions[stateseq]
            norm = sum(counts.values())
            for nextnucl in sorted(counts):
                count = counts[nextnucl]
                prob = count / norm
                yield stateseq, nextnucl, prob

    def __str__(self):
        output = 'State,Next,Prob'
        for initseq, prob in self.initial_probs:
            output += '\n{:s},,{:.4f}'.format(initseq, prob)
        for stateseq, nextnucl, prob in self.transition_probs:
            output += '\n{:s},{:s},{:.4f}'.format(stateseq, nextnucl, prob)
        return output


# -----------------------------------------------------------------------------
# Subcommand main methods
# -----------------------------------------------------------------------------

def train(args):
    model = NuclCompModel(order=args.order)
    for defline, sequence in parse_fasta(chain(*args.fasta)):
        model.train(sequence)
    print(model, file=args.out)


# -----------------------------------------------------------------------------
# Command-line interface
# -----------------------------------------------------------------------------

def train_parser(subparsers):
    subparser = subparsers.add_parser('train')
    subparser.add_argument('-o', '--out', type=argparse.FileType('w'),
                           metavar='FILE', help='write output to FILE')
    subparser.add_argument('-r', '--order', type=int, default=1, metavar='N',
                           help='order of the Markov model (nucleotide length '
                           'of the previous state); default is 1')
    subparser.add_argument('fasta', type=argparse.FileType('r'), nargs='+',
                           help='training sequence file(s) in Fasta format')


def get_parser():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='cmd', metavar='cmd')
    train_parser(subparsers)
    return parser


# -----------------------------------------------------------------------------
# Tests
# -----------------------------------------------------------------------------

def test_basic():
    model = NuclCompModel(order=1)
    assert model.order == 1

    with pytest.raises(AssertionError) as ae:
        model = NuclCompModel(order=0)

    with pytest.raises(AssertionError) as ae:
        model = NuclCompModel(order=-1)

    model = NuclCompModel(order=3)
    assert str(model) == 'State,Next,Prob'


def test_train():
    model = NuclCompModel(order=2)
    model.train('AAATAAATAAAT')
    assert sorted(model._transitions.keys()) == ['AA', 'AT', 'TA']
    assert list(model._initial_states.keys()) == ['AA']


def test_init_probs():
    model = NuclCompModel(order=2)
    model.train('AAATAAATAAAT')
    state, prob = next(model.initial_probs)
    assert state == 'AA'
    assert prob == 1.0


def test_trans_probs():
    model = NuclCompModel(order=2)
    model.train('AAATAAATAAAT')
    for state, nextnucl, prob in model.transition_probs:
        if state == 'AA' and nextnucl in ['A', 'T']:
            assert prob == 0.5, (state, nextnucl, prob)


def test_load():
    model = NuclCompModel(order=3)
    model_str = (
        'State,Next,Prob\n'
        'AAA,,1.0000\n'
        'AAA,A,0.0000\n'
        'AAA,C,0.0000\n'
        'AAA,G,0.0000\n'
        'AAA,T,1.0000\n'
        'AAT,A,1.0000\n'
        'AAT,C,0.0000\n'
        'AAT,G,0.0000\n'
        'AAT,T,0.0000\n'
        'ATA,A,1.0000\n'
        'ATA,C,0.0000\n'
        'ATA,G,0.0000\n'
        'ATA,T,0.0000\n'
        'TAA,A,1.0000\n'
        'TAA,C,0.0000\n'
        'TAA,G,0.0000\n'
        'TAA,T,0.0000\n'
    )
    model.load(iter(model_str.split('\n')))
    assert sorted(model._initial_states.keys()) == ['AAA']
    assert sorted(model._transitions.keys()) == ['AAA', 'AAT', 'ATA', 'TAA']


def test_load_headers():
    model = NuclCompModel(order=5)
    model_str = (
        'State,Next,Prob\n'
        'State,Next,Prob\n'
        'AAAAA,A,0.1000\n'
        'AAAAA,C,0.2000\n'
        'AAAAA,G,0.3000\n'
        'AAAAA,T,0.4000\n'
    )
    model.load(iter(model_str.split('\n')), num_header_lines=2)
    assert sorted(model._transitions.keys()) == ['AAAAA']

    model_str = (
        'AAAAA,A,0.1000\n'
        'AAAAA,C,0.2000\n'
        'AAAAA,G,0.3000\n'
        'AAAAA,T,0.4000\n'
    )
    model.load(iter(model_str.split('\n')), num_header_lines=0)
    assert sorted(model._transitions.keys()) == ['AAAAA']


def test_order():
    model = NuclCompModel(order=3)
    model_str = (
        'State,Next,Prob\n'
        'AA,A,0.25\n'
        'AA,C,0.25\n'
        'AA,G,0.25\n'
        'AA,T,0.25\n'
    )
    with pytest.raises(AssertionError) as ae:
        model.load(iter(model_str.split('\n')))


# -----------------------------------------------------------------------------
# OS-level CLI, not invoked when the module is imported
# -----------------------------------------------------------------------------

if __name__ == '__main__':
    mains = {
        'train': train,
    }

    args = get_parser().parse_args()
    mainmethod = mains[args.cmd]
    mainmethod(args)
