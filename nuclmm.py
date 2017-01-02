#!/usr/bin/env python
from __future__ import print_function, division
from collections import deque, defaultdict
from itertools import chain
import argparse
import re
import sys


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
    _initial_states = defaultdict(int)
    _transitions = defaultdict(nucldict)
    _nuclstate = deque()
    _order = int()

    def __init__(self, order=1):
        assert order > 0
        self._order = order

    def train(self, sequence):
        """
        Train the model on the given sequence.

        :param sequence: string-like object containing a nucleotide sequence
        """

        self._record_initial_state(sequence)
        for subseq in re.split('[^ACGT]+', sequence):
            for nucl, nextnucl in nuclpairs(subseq):
                self._nuclstate.append(nucl)
                if len(self._nuclstate) < self._order:
                    continue
                state = ''.join([n for n in self._nuclstate])
                self._transitions[state][nextnucl] += 1
                self._nuclstate.popleft()
            self._nuclstate = deque()  # Reset for next sequence/subsequence

    def _record_initial_state(self, sequence):
        init = sequence[:self._order]
        if re.search('[^ACGT]', init):
            message = ('Initial state "{}" includes an ambiguous/erroneous'
                       'nucleotide, ignoring'.format(init))
            print(message, file=sys.stderr)
        else:
            self._initial_states[init] += 1

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
# OS-level CLI, not invoked when the module is imported
# -----------------------------------------------------------------------------

if __name__ == '__main__':
    mains = {
        'train': train,
    }

    args = get_parser().parse_args()
    mainmethod = mains[args.cmd]
    mainmethod(args)
