#!/usr/bin/env python3

"""
Setup script for graphkernels package. Uses SWIG.
"""

from distutils import ccompiler
import pathlib

import numpy as np
import pkgconfig
import setuptools

import graphkernels

THIS_PATH = pathlib.Path('.')
GK_PATH = THIS_PATH / 'graphkernels'
CPP_PATH = GK_PATH / 'cppkernels'


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

SWIG_OPTS = ['-builtin', '-c++', '-O', '-py3', '-Wall']

try:
    _INFO = pkgconfig.parse('eigen3 python3')
except pkgconfig.PackageNotFoundError as e:
    e.message += """
        Missing `eigen3` library. Please install it using the
        package manager of your operating system.
        """
    raise

INCLUDE_DIRS = [np.get_include(), *_INFO['include_dirs'], str(CPP_PATH)]

LIBRARIES = _INFO['libraries']

CPP_SOURCES = [str(s) for s in CPP_PATH.glob('*.cpp')]


def main():
    _include_dir_flags = ccompiler.gen_preprocess_options(
        macros=[], include_dirs=INCLUDE_DIRS
    )
    setuptools.setup(
        ext_modules=[
            setuptools.Extension(
                '_graphkernels',
                sources=[str(GK_PATH / 'graphkernels.i'), *CPP_SOURCES],
                swig_opts=[*SWIG_OPTS, *_include_dir_flags],
                extra_compile_args=CPP_FLAGS,
                extra_link_args=CPP_FLAGS,
                include_dirs=INCLUDE_DIRS,
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
