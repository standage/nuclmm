# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of nuclmm (http://github.com/standage/nuclmm) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------
from __future__ import print_function
from itertools import chain
import nuclmm


def main(args):
    model = nuclmm.MarkovChain(order=args.order)
    for defline, sequence in nuclmm.parse_fasta(chain(*args.fasta)):
        model.train(sequence)
    print(model, file=args.out)
