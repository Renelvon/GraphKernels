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
CPP_DIR = path.join(GK_DIR, 'cppkernels')


def get_eigen_include_dir():
    try:
        cflags = pkgconfig.cflags('eigen3')
    except pkgconfig.PackageNotFoundError:
        print(
            """
            Missing `eigen3` library. Please install it using the
            package manager of your operating system.
            """,
            file=sys.stderr,
        )
        raise

    # Throw away the `-I` part; it is not part of the include directory.
    return cflags[2:]


CPP_FLAGS = [
    '-flto',
    '-g0',
    '-march=native',
    '-O3',
    '-std=c++14',
    '-Wall',
    '-Wextra',
    '-Wpedantic',
    '-Wno-sign-compare',
    '-Wno-unused-parameter',
    '-Wno-unused-variable',
]

SWIG_OPTS = ('-builtin', '-c++', '-O', '-py3', '-Wall')

LIBRARY_DIRS = [get_eigen_include_dir(), np.get_include(), CPP_DIR]

INCLUDE_DIR_FLAGS = tuple('-I{}'.format(l) for l in LIBRARY_DIRS)


def main():
    setuptools.setup(
        ext_modules=[
            setuptools.Extension(
                '_graphkernels',
                sources=[
                    path.join(GK_DIR, 'graphkernels.i'),  # Interface
                    path.join(CPP_DIR, 'connected_graphlet.cpp'),
                    path.join(CPP_DIR, 'graphlet.cpp'),
                    path.join(CPP_DIR, 'rest.cpp'),
                    path.join(CPP_DIR, 'wl.cpp'),
                ],
                swig_opts=(*SWIG_OPTS, *INCLUDE_DIR_FLAGS),
                extra_compile_args=CPP_FLAGS,
                extra_link_args=CPP_FLAGS,
                include_dirs=LIBRARY_DIRS,
                language='c++',
                optional=False,
            )
        ],
        version=graphkernels.__version__,
    )
    # Rest of options are specified in `setup.cfg`


if __name__ == '__main__':
    main()
