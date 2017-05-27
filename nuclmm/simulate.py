# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of nuclmm (http://github.com/standage/nuclmm) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function
import argparse
import random
import nuclmm


def subparser(subparsers):
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
    subparser.add_argument('-i', '--ignore-inits', action='store_true',
                           help='ignore initial state probabilities from the '
                           'given model and generate initial states randomly')
    subparser.add_argument('modelfile', type=argparse.FileType('r'))


def main(args):
    if args.seed:
        random.seed(args.seed)
    model = nuclmm.MarkovChain(order=args.order)
    model.load(args.modelfile)
    simulator = model.simulate(args.numseqs, args.seqlen, args.ignore_inits)
    for i, sequence in enumerate(simulator, 1):
        print('>seq', i, '\n', sequence, sep='', file=args.out)
