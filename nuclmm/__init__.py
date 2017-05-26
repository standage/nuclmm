# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of nuclmm (http://github.com/standage/nuclmm) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

try:
    import __builtin__ as builtins
except:  # pragma: no cover
    import builtins
import argparse
import pkg_resources
from nuclmm import train
from nuclmm import simulate
from nuclmm.markovchain import MarkovChain


from ._version import get_versions
__version__ = get_versions()['version']
del get_versions


def open(filename, mode):
    if mode not in ['r', 'w']:
        raise ValueError('invalid mode "{}"'.format(mode))
    if filename in ['-', None]:
        filehandle = sys.stdin if mode == 'r' else sys.stdout
        return filehandle
    openfunc = builtins.open
    if filename.endswith('.gz'):
        openfunc = gzopen
        mode += 't'
    return openfunc(filename, mode)


def data_file(basename):
    datadir = pkg_resources.resource_filename('nuclmm', 'tests/data')
    return datadir + '/' + basename


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


def cli_parser():
    desc = ('Use arbitrary-order Markov chains to model nucleotide composition'
            ' and simulate random sequences.')
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-v', '--version', action='version',
                        version='nuclmm v{}'.format(__version__))
    parser._positionals.title = 'Subcommands'
    parser._optionals.title = 'Global arguments'

    subcommandstr = '", "'.join(['train, simulate'])
    subparsers = parser.add_subparsers(dest='cmd', metavar='cmd',
                                       help=subcommandstr)
    train.subparser(subparsers)
    simulate.subparser(subparsers)
    return parser
