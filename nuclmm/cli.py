# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of nuclmm (http://github.com/standage/nuclmm) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

import argparse
import nuclmm


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
    desc = ('Use arbitrary-order Markov chains to model nucleotide composition'
            ' and simulate random sequences.')
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-v', '--version', action='version',
                        version='nuclmm v{}'.format(nuclmm.__version__))
    parser._positionals.title = 'Subcommands'
    parser._optionals.title = 'Global arguments'

    subcommandstr = '", "'.join(['train, simulate'])
    subparsers = parser.add_subparsers(dest='cmd', metavar='cmd',
                                       help=subcommandstr)
    train_parser(subparsers)
    simulate_parser(subparsers)
    return parser
