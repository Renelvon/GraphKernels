from setuptools import setup

import graphkernels


def main():
    setup(
        author='Elisabetta Ghisu',
        author_email='elisabetta.ghisu@bsse.ethz.ch',
        classifiers=[
            'Development Status :: 4 - Beta',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3',
            'Topic :: Scientific/Engineering :: Mathematics',
        ],
        description='Package for computing graph kernels',
        install_requires=['numpy', 'python-igraph', 'GKextCPy'],
        license='ETH Zurich',
        long_description='', # TODO: Fill me!
        name='graphkernels',
        packages=['graphkernels'],
        package_data={'graphkernels': ['data.mutag']},
        python_requires=">=2.7, !=3.0.*, !=3.1.*, !=3.2.*, !=3.3.*",
        version=graphkernels.__version__,
        url='https://github.com/BorgwardtLab/GraphKernels',
        zip_safe=False
    )


if __name__ == '__main__':
    main()
