# graphkernels

[![Build Status](https://travis-ci.org/Renelvon/GraphKernels.svg?branch=master)](https://travis-ci.org/Renelvon/GraphKernels)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

`graphkernels` is a Python package for computing various graph kernels. The
[backend](https://github.com/mahito-sugiyama/graph-kernels) is coded in C++ and
exposed to Python through [SWIG](http://www.swig.org/).

The Python and R packages are described at:

- M. Sugiyama, M.E. Ghisu, F. Llinares-López and K. Borgwardt. **graphkernels:
  R and Python packages for graph comparison**. Bioinformatics, 2017.

The paper can be found
[here](https://academic.oup.com/bioinformatics/article/34/3/530/4209994/)

## Installation

### Installing dependencies

`graphkernels` requires the following dependencies to be pre-installed in your
local environment:

- Python bindings for [igraph](https://pypi.org/project/python-igraph/)
- [NumPy](https://pypi.org/project/numpy/)
- [Eigen3](http://eigen.tuxfamily.org/)

During setup, you will also need the following dependencies:

- a C++ compiler (e.g. [gcc](http://gcc.gnu.org),
  [XCode](https://developer.apple.com/xcode/))
- [pkg-config](https://www.freedesktop.org/wiki/Software/pkg-config/) and the
  respective [Python bindings](https://pypi.org/project/pkgconfig/)
- [SWIG](http://www.swig.org)

We recommend the following steps for installing the dependencies (except for
the C++ compiler):

- On Ubuntu 18.04 LTS:

    `$ sudo apt-get install libeigen3-dev python3-igraph python3-numpy
    python3-pkgconfig swig`

Some of the packages might not be available in older distributions.

- On MacOSX:

    `$ brew install eigen pkg-config`

Some dependencies may be skipped without problem, because they are
automatically downloaded and installed by `pip`.

To install the package through [PyPI](https://pypi.org/), type:

    $ pip3 install graphkernels

It is recommended that version `10` or later of `pip` is used. To see which
version is installed, type:

    $ pip3 --version

If you are using an older version, you can locally upgrade `pip` by typing:

    $ python3 -m pip install --user --upgrade-pip

and retrying the `graphkernels` installation.

Alternatively, the package can be build from source. After downloading the
source code from GitHub

    $ git clone https://github.com/Renelvon/GraphKernels.git

users can use the `setup.py` script to install the package

    $ cd GraphKernels
    $ python3 setup.py build_ext
    $ python3 setup.py install --user

You should also make sure that you're installing the latest release of our
package, in case you've had a previous version installed. To make sure the
extension and package are not taken from your cache, you can use the
`--no-cache-dir` option and install the package as:

    $ pip3 --no-cache-dir install graphkernels

## Usage

The `tutorial` folder should help with getting started using the package.
There, you can also find an example script for computing graph kernels through
our package on a benchmark dataset.

## Compatibility

This version of `graphkernels` aims to be fast, user-friendly and
backwards-compatible (in that order). Accordingly, support for Python 2 in
`graphkernels` has been dropped because it is reaching
[End-Of-Life](https://pythonclock.org/) very soon. If you care about speed
enough to use a C++ backend for your project, you should care as much about
future portability and give serious thought before continuing development on
Python 2. Most of the [ML projects](https://python3statement.org/) will drop
support for Python 2 soon.

## Citation

If you use the `graphkernels` package in your projects please cite our work as
follows:

``` @article{Sugiyama-2017-Bioinformatics,
author = {Sugiyama, Mahito and Ghisu, M. Elisabetta and Llinares-López, Felipe and Borgwardt, Karsten},
title = {graphkernels: R and Python packages for graph comparison},
journal = {Bioinformatics},
volume = {34},
number = {3},
pages = {530--532},
year = {2017},
doi = {10.1093/bioinformatics/btx602},
URL = {http://dx.doi.org/10.1093/bioinformatics/btx602},
} ```

If our project has made your life easier, we will be happy to receive your
feedback, through e-mail or an issue on GitHub.
