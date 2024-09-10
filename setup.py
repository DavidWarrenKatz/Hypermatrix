#!/usr/bin/env python3

from setuptools import setup, find_packages

setup(
    name="Hypermatrix",
    version="0.1.0",
    packages=find_packages(include=['hypermatrix', 'hypermatrix.*']),
    package_data={
        'hypermatrix': [
            'utilities/**/*',  # file i want to include
        ],
    },
    author="David Warren Katz",
    author_email="davidkatz02@gmail.com",
    description=(
        "Hypermatrix is a command-line tool for integrating multi-omics data and "
        "analyzing epigenetic data using advanced tensor techniques. It features commands "
        "like 'ABcluster' for single-cell clustering and A/B compartment analysis, and "
        "'differentiate_chromosomes' for distinct A/B compartment calls across homologous chromosomes. "
        "Hypermatrix leverages Non-Negative Tensor Factorization (NTF) to preserve critical interaction "
        "information, offering a powerful solution for complex multi-omics integration and epigenetic analysis."
    ),
    url="https://github.com/DavidWarrenKatz/Hypermatrix",
    entry_points={
        'console_scripts': [
            'hypermatrix=hypermatrix.main:main',
        ],
    },
    install_requires=[
        'numpy',
        'matplotlib',
        'seaborn',
        'requests',
        'h5py',
        'scipy',
        'pyBigWig',
        'cooler',
        'cytoolz',
        'future',
        'intervaltree',
        'msgpack',
        'pytest',
        'scikit-image',
        'scikit-learn',
        'tables',
    ],
    python_requires='>=3.8',
)

