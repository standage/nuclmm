# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of nuclmm (http://github.com/standage/nuclmm) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------
import nuclmm


def main():
    mains = {
        'train': nuclmm.train.main,
        'simulate': nuclmm.simulate.main,
    }

    args = nuclmm.cli.get_parser().parse_args()
    mainmethod = mains[args.cmd]
    mainmethod(args)
