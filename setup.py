#!/usr/bin/env python3

"""
Setup script for graphkernels package. Uses SWIG.
"""

from os import path
import sys

import numpy as np
import pkgconfig
import setuptools

import graphkernels

THIS_DIR = path.dirname(__file__)
GK_DIR = path.join(THIS_DIR, 'graphkernels')


def get_eigen_include_dir():
    try:
        cflags = pkgconfig.cflags('eigen3')
    except pkgconfig.PackageNotFoundError:
        print(
            """
            Missing `eigen3` library. Please install it using the
            package manager of your operating system.
            """,
            file=sys.stderr
        )
        raise

    # Throw away the `-I` part; it is not part of the include directory.
    return cflags[2:]


def main():
    cppflags = [
        '-flto',
        '-g0',
        '-march=native',
        '-O3',
        '-std=c++14',
        '-Wall',
        '-Wextra',
        '-Wpedantic'
    ]

    setuptools.setup(
        ext_modules=[
            setuptools.Extension(
                '_graphkernels',
                sources=[
                    # Interface file
                    path.join(GK_DIR, 'graphkernels.i'),

                    # Implementation file
                    path.join(GK_DIR, 'graphkernels.cpp'),
                ],
                swig_opts=['-c++', '-Wall', '-builtin', '-O', '-py3'],
                extra_compile_args=cppflags,
                extra_link_args=cppflags,
                include_dirs=[
                    get_eigen_include_dir(),
                    np.get_include(),
                ]
            )
        ],
        # NOTE: The following option may produce a harmless warning when
        # building package using setuptools < 40.*
        long_description_content_type='text/markdown',
        version=graphkernels.__version__,
    )
    # Rest of options are specified in `setup.cfg`


if __name__ == '__main__':
    main()
