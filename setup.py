#!/usr/bin/env python3

"""
Setup script for graphkernels package. Uses SWIG.
"""

from os import path

import numpy as np
import pkgconfig
import setuptools

import graphkernels

THIS_DIR = path.dirname(__file__)
GK_DIR = path.join(THIS_DIR, 'graphkernels')
CPP_DIR = path.join(GK_DIR, 'cppkernels')


CPP_FLAGS = [
    '-flto',
    '-g0',
    '-march=native',
    '-O3',
    '-std=c++14',
    '-Wall',
    '-Werror',
    '-Wextra',
    '-Wl,-z,defs',
    '-Wpedantic',
    '-Wno-sign-compare',
    '-Wno-unused-parameter',
    '-Wno-unused-variable',
]

SWIG_OPTS = ('-builtin', '-c++', '-O', '-py3', '-Wall')

try:
    _INFO = pkgconfig.parse('eigen3 python3')
except pkgconfig.PackageNotFoundError as e:
    e.message += """
        Missing `eigen3` library. Please install it using the
        package manager of your operating system.
        """
    raise

LIBRARY_DIRS = [np.get_include(), *_INFO['include_dirs'], CPP_DIR]

LIBRARIES = _INFO['libraries']

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
                libraries=LIBRARIES,
                language='c++',
                optional=False,
            )
        ],
        version=graphkernels.__version__,
    )
    # Rest of options are specified in `setup.cfg`


if __name__ == '__main__':
    main()
