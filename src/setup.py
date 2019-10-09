#!/usr/bin/env python

"""
setup.py file for SWIG
"""

from os import path
from setuptools import setup, Extension

import pkgconfig
import numpy


def get_include_dirs():

    if not pkgconfig.exists('eigen3'):
        raise Exception('Missing `eigen3` library. Please install it using the package manager of your operating system')

    np_include_dir = numpy.get_include()

    # We have to throw away the `-I` part of the `pkgconfig` output
    # because it is not part of the include directory.
    eigen3_include_dir = pkgconfig.cflags('eigen3')[2:]

    return [np_include_dir, eigen3_include_dir]


here = path.abspath(path.dirname(__file__))
gk_dir = path.join(here, 'GKextCPy')


GKextCPy_module = Extension(
    '_GKextCPy',
    sources=[
        path.join(gk_dir, 'GKextCPy_wrap.cxx'),
        path.join(gk_dir, 'GKextCPy.cpp'),
    ],
    swig_opts=['-c++'],
    extra_compile_args=['-std=c++11', '-O3'],
    include_dirs=get_include_dirs()
)


setup(
    author='Elisabetta Ghisu',
    author_email='elisabetta.ghisu@bsse.ethz.ch',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Programming Language :: C++',
        'Programming Language :: Python :: 3 :: Only',
        'Topic :: Scientific/Engineering :: Mathematics',
    ],
    description="""
        Graph Kernels: building the extension Python module.
        This is a wrapper package from C++ to Python.
    """,
    ext_modules=[GKextCPy_module],
    install_requires=['pkgconfig', 'numpy'],
    license='ETH Zurich',
    long_description='', # TODO: Fill me!
    name='GKextCPy',
    package_dir={'': 'GKextCPy'},
    py_modules=['GKextCPy'],
    python_requires='>=3.4',
    setup_requires=['pkgconfig', 'numpy'],
    url='https://github.com/BorgwardtLab/GraphKernels',
    version='0.5.0', # TODO: Import version from init file.
    zip_safe=False
)
