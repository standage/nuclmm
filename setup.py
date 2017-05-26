#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (c) 2017 The Regents of the University of California
#
# This file is part of nuclmm (http://github.com/standage/nuclmm) and is
# licensed under the MIT license: see LICENSE.
# -----------------------------------------------------------------------------

from setuptools import setup
import glob
import versioneer


desc = ('Use arbitrary-order Markov chains to model nucleotide composition '
        'and simulate random sequences')
setup(name='nuclmm',
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
      description=desc,
      url='https://github.com/standage/nuclmm',
      author='Daniel Standage',
      author_email='daniel.standage@gmail.com',
      license='MIT',
      packages=['nuclmm', 'nuclmm.tests'],
      package_data={
          'nuclmm': ['nuclmm/tests/data/*']
      },
      entry_points={'console_scripts': ['nuclmm = nuclmm.__main__:main']},
      classifiers=[
          'Development Status :: 4 - Beta',
          'Environment :: Console',
          'License :: OSI Approved :: BSD License',
          'Programming Language :: Python :: 2.7',
          'Programming Language :: Python :: 3.4',
          'Programming Language :: Python :: 3.5',
          'Topic :: Scientific/Engineering :: Bio-Informatics'
      ],
      zip_safe=True)
