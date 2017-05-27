# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of nuclmm (http://github.com/standage/nuclmm) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from __future__ import print_function
import nuclmm
import pytest


def test_help(capsys):
    with pytest.raises(SystemExit):
        nuclmm.cli_parser().parse_args(['-h'])
    out, err = capsys.readouterr()
    assert 'show this help message and exit' in out


def test_version(capsys):
    with pytest.raises(SystemExit):
        nuclmm.cli_parser().parse_args(['-v'])
    out, err = capsys.readouterr()
    assert nuclmm.__version__ in out or nuclmm.__version__ in err
