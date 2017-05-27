# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of nuclmm (http://github.com/standage/nuclmm) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------
import nuclmm


def main(args=None):
    mains = {
        'train': nuclmm.train.main,
        'simulate': nuclmm.simulate.main,
    }

    if args is None:  # pragma: no cover
        args = nuclmm.cli_parser().parse_args()

    if args.cmd is None:  # pragma: no cover
        nuclmm.cli_parser().parse_args(['-h'])

    args = nuclmm.cli_parser().parse_args()
    mainmethod = mains[args.cmd]
    mainmethod(args)
