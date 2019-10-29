# graphkernels

[![Build Status](https://travis-ci.org/Renelvon/GraphKernels.svg?branch=master)](https://travis-ci.org/Renelvon/GraphKernels)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![Code style: Black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

`graphkernels` is a Python package for computing various graph kernels. The
[backend](https://github.com/mahito-sugiyama/graph-kernels) is coded in C++ and
exposed to Python through [SWIG](http://www.swig.org/).

The Python and R packages are described at:

- M. Sugiyama, M.E. Ghisu, F. Llinares-López and K. Borgwardt. **graphkernels:
  R and Python packages for graph comparison**. Bioinformatics, 2017.

The paper can be found
[here](https://academic.oup.com/bioinformatics/article/34/3/530/4209994/).

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

```sh
$ sudo apt-get install libeigen3-dev python3-igraph python3-numpy python3-pkgconfig swig
```

Some of the packages might not be available in older distributions.

- On MacOSX:

```sh
$ brew install eigen pkg-config
```

### Installing `graphkernels` from PyPI

To install the package through [PyPI](https://pypi.org/), type:

```sh
$ pip3 install graphkernels --user
```

It is recommended that version `10` or later of `pip` is used; newer versions
of `pip` can handle missing dependencies better, by automatically downloading
and installing them together with `graphkernels` (this won't work for the Eigen
library). To see which version of `pip` is installed, type:

```sh
$ pip3 --version
```

If you are using an older version, you can locally upgrade `pip` by typing:

```
$ python3 -m pip install --user --upgrade-pip
```

before attempting the `graphkernels` installation.

In case you already have an older version of `graphkernels` installed, you can
use the `--no-cache-dir` option so that your local package cache is ignored and
the lastest package is downloaded:

```sh
$ pip3 --no-cache-dir install graphkernels --user
```

### Build-Installing `graphkernels` from main repository
Alternatively, the package can be built from source as follows:

1. Download the latest source code from GitHub (please *only* use the `master`
   branch when doing so, other branches are not guaranteed to contain a correct
   build):

```sh
$ git clone https://github.com/Renelvon/GraphKernels.git
```

2. Build the package; this will also compile the C++ backend and generate the
   SWIG wrapper:

```sh
$ cd GraphKernels
$ python3 setup.py build
```

3. Install the package for this user only (showed by `--user`; installing the
   package system-wide might require `sudo -H`):

```sh
$ python3 setup.py install --user
```

## Usage

The `Tutorial` folder should help with getting started using the package.
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

```
@article{Sugiyama-2017-Bioinformatics,
author = {Sugiyama, Mahito and Ghisu, M. Elisabetta and Llinares-López, Felipe and Borgwardt, Karsten},
title = {graphkernels: R and Python packages for graph comparison},
journal = {Bioinformatics},
volume = {34},
number = {3},
pages = {530--532},
year = {2017},
doi = {10.1093/bioinformatics/btx602},
URL = {http://dx.doi.org/10.1093/bioinformatics/btx602},
}
```

If our project has made your life easier, we will be happy to receive your
feedback, through e-mail or an issue on GitHub.
