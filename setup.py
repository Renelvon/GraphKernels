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

README = 'README.md'
THIS_DIR = path.dirname(__file__)
GK_DIR = path.join(THIS_DIR, 'graphkernels')
README_PATH = path.join(THIS_DIR, README)


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
    with open(README_PATH, 'r') as f_readme:
        long_description = f_readme.read()

    setuptools.setup(
        author='Elisabetta Ghisu',
        author_email='elisabetta.ghisu@bsse.ethz.ch',
        classifiers=[
            'Development Status :: 4 - Beta',
            'Programming Language :: C++',
            'Programming Language :: Python :: 3 :: Only',
            'Programming Language :: Python :: 3.5',
            'Programming Language :: Python :: 3.6',
            'Programming Language :: Python :: 3.7',
            'Programming Language :: Python :: 3.8',
            'Topic :: Scientific/Engineering :: Mathematics',
        ],
        description='Package for computing graph kernels',
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
                extra_compile_args=['-std=c++11', '-O3'],
                include_dirs=[
                    get_eigen_include_dir(),
                    np.get_include(),
                ]
            )
        ],
        include_package_data=True,
        install_requires=[
            'numpy>=1.11',
            'pkgconfig',
            'python-igraph'
        ],
        license='ETH Zurich',
        long_description=long_description,
        # NOTE: The following option may produce a harmless warning when
        # building package using setuptools < 40.*
        long_description_content_type='text/markdown',
        name='graphkernels',
        packages=['graphkernels'],
        python_requires='>=3.4',
        setup_requires=[
            'numpy>=1.11',
            'pkgconfig'
        ],
        tests_require=['cpplint', 'pylint'],
        url='https://github.com/BorgwardtLab/GraphKernels',
        version=graphkernels.__version__,
        zip_safe=False
    )


if __name__ == '__main__':
    main()
