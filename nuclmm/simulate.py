# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of nuclmm (http://github.com/standage/nuclmm) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------
import nuclmm


def main(args):
    if args.seed:
        random.seed(args.seed)
    model = nuclmm.MarkovChain(order=args.order)
    model.load(args.modelfile)
    for i, sequence in enumerate(model.simulate(args.numseqs, args.seqlen)):
        print('>seq', i, '\n', sequence, sep='', file=args.out)
