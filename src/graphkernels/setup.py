from setuptools import setup

import graphkernels


def main():
    setup(
        name='graphkernels',
        version=graphkernels.__version__,
        description='Package for computing graph kernels',
        url='https://github.com/BorgwardtLab/GraphKernels',
        author='Elisabetta Ghisu',
        author_email='elisabetta.ghisu@bsse.ethz.ch',
        license='ETH Zurich',
        packages=['graphkernels'],
        install_requires=['GKextCPy', 'python-igraph', 'numpy'],
        package_data={'graphkernels': ['data.mutag']},
    )


if __name__ == '__main__':
    main()
