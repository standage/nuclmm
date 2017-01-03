#!/usr/bin/env python
from __future__ import print_function, division
from collections import deque, defaultdict
from itertools import chain
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO
import argparse
import random
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


def draw(dist):
    rn = random.random()
    test = 0.0
    for key in sorted(dist.keys()):
        test += dist[key]
        if rn < test:
            return key
    return key


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
        self._normed = False
        self._initial_states = defaultdict(int)
        self._transitions = defaultdict(nucldict)
        self._nuclstate = deque()
        self._order = int(order)

    def train(self, sequence):
        """
        Train the model on the given sequence.

        :param sequence: string-like object containing a nucleotide sequence
        """
        self._normed = False
        self._record_initial_state(sequence)
        for subseq in re.split('[^ACGT]+', sequence.upper()):
            for nucl, nextnucl in nuclpairs(subseq):
                self._nuclstate.append(nucl)
                if len(self._nuclstate) < self._order:
                    continue
                self._transitions[self.state][nextnucl] += 1
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
        self._normed = True

    def simulate(self, numseqs, seqlen, ignore_inits=False):
        for _ in range(numseqs):
            seq = deque()
            self._nuclstate = deque()
            initseq = draw(self._initial_states)
            for n in initseq:
                self._nuclstate.append(n)
                seq.append(n)
            while len(seq) < seqlen:
                nextnucl = draw(self._transitions[self.state])
                seq.append(nextnucl)
                self._nuclstate.append(nextnucl)
                self._nuclstate.popleft()
            yield ''.join([n for n in seq])

    def normalize(self):
        if self._normed:
            return

        probs = dict()
        norm = sum(self._initial_states.values())
        for initseq in self._initial_states:
            count = self._initial_states[initseq]
            prob = count / norm
            probs[initseq] = prob
        self._initial_states = probs

        probs = defaultdict(nucldict)
        for stateseq in self._transitions:
            counts = self._transitions[stateseq]
            norm = sum(counts.values())
            for nextnucl in sorted(counts):
                count = counts[nextnucl]
                prob = count / norm
                probs[stateseq][nextnucl] = prob
        self._transitions = probs

        self._normed = True

    def _record_initial_state(self, sequence):
        init = sequence[:self._order]
        if re.search('[^ACGT]', init):
            message = ('Initial state "{}" includes an ambiguous/erroneous'
                       'nucleotide, ignoring'.format(init))
            print(message, file=sys.stderr)
        else:
            self._initial_states[init] += 1

    @property
    def state(self):
        return ''.join([n for n in self._nuclstate])

    @property
    def order(self):
        return self._order

    @property
    def initial_probs(self):
        self.normalize()
        for initseq in sorted(self._initial_states):
            prob = self._initial_states[initseq]
            yield initseq, prob

    @property
    def transition_probs(self):
        self.normalize()
        for stateseq in sorted(self._transitions):
            probs = self._transitions[stateseq]
            for nextnucl in sorted(probs):
                prob = probs[nextnucl]
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


def simulate(args):
    if args.seed:
        random.seed(args.seed)
    model = NuclCompModel(order=args.order)
    model.load(args.modelfile)
    for i, sequence in enumerate(model.simulate(args.numseqs, args.seqlen)):
        print('>seq', i, '\n', sequence, sep='', file=args.out)


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


def simulate_parser(subparsers):
    subparser = subparsers.add_parser('simulate')
    subparser.add_argument('-o', '--out', type=argparse.FileType('w'),
                           metavar='FILE', help='write output to FILE')
    subparser.add_argument('-r', '--order', type=int, default=1, metavar='N',
                           help='order of the Markov model (nucleotide length '
                           'of the previous state); default is 1')
    subparser.add_argument('-n', '--numseqs', type=int, metavar='N',
                           default=1000, help='number of sequences to produce;'
                           ' default is 1000')
    subparser.add_argument('-l', '--seqlen', type=int, metavar='L',
                           default=100, help='sequence length; default is 100')
    subparser.add_argument('-s', '--seed', type=int, metavar='S',
                           help='seed for random number generator')
    subparser.add_argument('modelfile', type=argparse.FileType('r'))


def get_parser():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='cmd', metavar='cmd')
    train_parser(subparsers)
    simulate_parser(subparsers)
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
    with open('aaat.fa', 'r') as infile:
        for defline, sequence in parse_fasta(infile):
            model.train(sequence)
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
    with open('aaat.mm', 'r') as infile:
        model.load(infile)
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


def test_order_mismatch():
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


def test_train_cli():
    with open('aaat.fa', 'r') as seqfile, open('aaat.mm', 'r') as modelfile:
        args = type('', (), {})()
        args.out = StringIO()
        args.order = 3
        args.fasta = [seqfile]
        train(args)

        model = modelfile.read().strip()
        modeltest = args.out.getvalue().strip()
        assert modeltest == model


def test_simulate_cli():
    with open('aaat.mm', 'r') as mfile, open('aaat-pass2.fa', 'r') as seqfile:
        args = type('', (), {})()
        args.out = StringIO()
        args.order = 3
        args.numseqs = 5
        args.seqlen = 50
        args.seed = 42
        args.modelfile = mfile
        simulate(args)

        seqs = seqfile.read().strip()
        seqstest = args.out.getvalue().strip()
        assert seqstest == seqs


# -----------------------------------------------------------------------------
# OS-level CLI, not invoked when the module is imported
# -----------------------------------------------------------------------------

if __name__ == '__main__':
    mains = {
        'train': train,
        'simulate': simulate,
    }

    args = get_parser().parse_args()
    mainmethod = mains[args.cmd]
    mainmethod(args)
