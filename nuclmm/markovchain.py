#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of nuclmm (http://github.com/standage/nuclmm) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

"""Nth order Markov chain as a model of nucleotide composition."""

from __future__ import print_function, division
from collections import deque, defaultdict
import json
import random
import re
import sys


def draw(dist):
    """
    Draw a key randomly from the distribution.

    The distribution is stored as a dictionary, where the keys are items of
    interest and the values are the associated probability of drawing a
    particular item.
    """
    rand = random.random()
    test = 0.0
    key = None
    for key in sorted(dist):
        test += dist[key]
        if rand < test:
            return key
    return key


def nuclpairs(sequence):
    """Generator yielding overlapping pairs of nucleotides from a sequence."""
    if len(sequence) < 2:
        return
    prev = sequence[0]
    for nucl in sequence[1:]:
        yield prev, nucl
        prev = nucl


def nucldict():
    """Initialize an empty dictionary of nucleotide counts or probabilities."""
    return {n: 0 for n in 'ACGT'}


class MarkovChainOrderMismatch(ValueError):
    pass


class MarkovChainNormalizationError(ValueError):
    pass


class MarkovChain(object):
    """
    Class for training models of nucleotide composition.

    In this context an Nth order Markov chain is defined as follows: first, a
    set of *transition probabilities*, i.e., the probability that the next
    nucleotide will be *X* given that the current state (*N* nucleotides) is
    *Y*; and second, a set of *initial state probabilities*, i.e., the
    probability of observing *Z* as the first *N* nucleotides of a given
    sequence.

    This class supports models of arbitrary order (N >= 1), although memory
    requirements grow exponentially with respect to N and will likely become
    intractable for N=10 or thereabouts.

    Probabilities are computed by observing invoking the **train** command on
    one or more training sequences.
    """

    def __init__(self, order=1):
        assert order > 0
        self._normed = False
        self._initial_states = defaultdict(int)
        self._transitions = defaultdict(nucldict)
        self._nuclstate = deque()
        self._order = int(order)

    def __str__(self):
        return json.dumps([self._initial_states, self._transitions], indent=4)

    def __eq__(self, other):
        if self.order != other.order:
            return False
        if set(list(self._transitions)) != set(list(other._transitions)):
            return False
        for state in self._transitions:
            if state not in other._transitions:
                return False
            for nucl in 'ACGT':
                selfprob = self._transitions[state][nucl]
                otherprob = other._transitions[state][nucl]
                if abs(selfprob - otherprob) > 0.0001:
                    return False
        return True

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

    def save(self, outstream):
        """Save the Markov chain to a file, in JSON format."""
        print(self, file=outstream)

    def load(self, instream):
        """Load the Markov chain from a JSON file."""
        data = json.load(instream)
        assert len(data) == 2
        self._initial_states = data[0]
        self._transitions = data[1]
        self._normed = True
        self.validate()

    def validate(self):
        allinits = 0.0
        for state, prob in self.initial_probs:
            if len(state) != self.order:
                message = 'state has length {:d} '.format(len(state))
                message += 'but model declared {:d}-order'.format(self.order)
                raise MarkovChainOrderMismatch(message)
            allinits += prob
        if abs(allinits - 1.0) > 0.0001:
            message = 'initial state probabilites should sum to 1.0'
            message += '; got {:.4f}'.format(allinits)
            raise MarkovChainNormalizationError(message)

        for state in self._transitions:
            probsum = 0.0
            for nextnucl in self._transitions[state]:
                probsum += self._transitions[state][nextnucl]
            if abs(probsum - 1.0) > 0.0001:
                message = 'state transition probabilities for {}'.format(state)
                message += ' should sum to 1.0; got {:.4f}'.format(probsum)
                raise MarkovChainNormalizationError(message)

    def random_init(self):
        """Generate a random initial state sequence."""
        return ''.join(random.choice('ACGT') for _ in range(self.order))

    def simulate(self, numseqs, seqlen, ignore_inits=False):
        """Simulate a sequence from the Markov chain."""
        self.normalize()
        for _ in range(numseqs):
            seq = deque()
            self._nuclstate = deque()

            if ignore_inits:
                initseq = self.random_init()
            else:
                initseq = draw(self._initial_states)
            for nucl in initseq:
                self._nuclstate.append(nucl)
                seq.append(nucl)

            while len(seq) < seqlen:
                nextnucl = draw(self._transitions[self.state])
                seq.append(nextnucl)
                self._nuclstate.append(nextnucl)
                self._nuclstate.popleft()
            yield ''.join([n for n in seq])

    def normalize(self):
        """
        Normalize counts to probabilities.

        Initial states and transition states are initially read in as counts.
        This function normalizes these counts to probabilities.
        """
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
            message = ('initial state "{}" includes an ambiguous/erroneous'
                       'nucleotide, ignoring'.format(init))
            print('[nuclmm::MarkovChain] Warning:', message, file=sys.stderr)
        else:
            self._initial_states[init] += 1

    @property
    def state(self):
        """String representation of the current state."""
        return ''.join([n for n in self._nuclstate])

    @property
    def order(self):
        """The order of the Markov chain."""
        return self._order

    @property
    def initial_probs(self):
        """Generator for initial state probabilities."""
        self.normalize()
        for initseq in sorted(self._initial_states):
            prob = self._initial_states[initseq]
            yield initseq, prob

    @property
    def transition_probs(self):
        """Generator for state transition probabilities."""
        self.normalize()
        for stateseq in sorted(self._transitions):
            probs = self._transitions[stateseq]
            for nextnucl in sorted(probs):
                prob = probs[nextnucl]
                yield stateseq, nextnucl, prob
