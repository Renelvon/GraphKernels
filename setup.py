#!/usr/bin/env python3

"""
Setup script for GKextCPy package. Uses SWIG.
"""

from os import path
import sys

import numpy as np
import pkgconfig
from setuptools import setup, Extension

import GKextCPy


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


GK_DIR = path.join(path.dirname(__file__), 'GKextCPy')

GKextCPy_module = Extension(
    '_GKextCPy',
    sources=[
        # Interface file
        path.join(GK_DIR, 'GKextCPy.i'),

        # Implementation file
        path.join(GK_DIR, 'GKextCPy.cpp'),
    ],
    swig_opts=['-c++', '-Wall', '-builtin', '-O', '-py3'],
    extra_compile_args=['-std=c++11', '-O3'],
    include_dirs=[
        get_eigen_include_dir(),
        np.get_include(),
    ]
)


def main():
    setup(
        author='Elisabetta Ghisu',
        author_email='elisabetta.ghisu@bsse.ethz.ch',
        classifiers=[
            'Development Status :: 4 - Beta',
            'Programming Language :: C++',
            'Programming Language :: Python :: 3 :: Only',
            'Topic :: Scientific/Engineering :: Mathematics',
        ],
        description='Package for computing graph kernels',
        ext_modules=[GKextCPy_module],
        include_package_data=True,
        install_requires=['numpy', 'pkgconfig', 'python-igraph'],
        license='ETH Zurich',
        long_description='', # TODO: Fill me!
        name='GKextCPy',
        packages=['GKextCPy'],
        python_requires='>=3.4',
        setup_requires=['pkgconfig', 'numpy'],
        tests_require=['cpplint', 'pylint'],
        url='https://github.com/BorgwardtLab/GraphKernels',
        version=GKextCPy.__version__,
        zip_safe=False
    )


if __name__ == '__main__':
    main()
