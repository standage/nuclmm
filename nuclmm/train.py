# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of nuclmm (http://github.com/standage/nuclmm) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------
from __future__ import print_function
import argparse
from itertools import chain
import nuclmm


def subparser(subparsers):
    subparser = subparsers.add_parser('train')
    subparser.add_argument('-o', '--out', type=argparse.FileType('w'),
                           metavar='FILE', help='write output to FILE')
    subparser.add_argument('-r', '--order', type=int, default=1, metavar='N',
                           help='order of the Markov model (nucleotide length '
                           'of the previous state); default is 1')
    subparser.add_argument('fasta', nargs='+',
                           help='training sequence file(s) in Fasta format')


def main(args):
    fastas = [nuclmm.open(f, 'r') for f in args.fasta]
    model = nuclmm.MarkovChain(order=args.order)
    for defline, sequence in nuclmm.parse_fasta(chain(*fastas)):
        model.train(sequence)
    model.normalize()
    print(model, file=args.out)
