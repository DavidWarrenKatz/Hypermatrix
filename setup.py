#!/usr/bin/env python3

from setuptools import setup, find_packages

setup(
    name="Hypermatrix",
    version="0.1.0",
    packages=find_packages(include=['hypermatrix', 'hypermatrix.*']),
    package_data={
        'hypermatrix': [
            'hypermatrix/utilities/**/*',  # Include all files in the utilities directory
        ],
    },
    include_package_data=True,
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
        'numpy>=1.21.0',
        'matplotlib>=3.4.0',
        'seaborn>=0.11.0',
        'requests>=2.25.0',
        'h5py>=3.3.0',
        'scipy>=1.7.0',
        'pyBigWig==0.3.22',  
    ],
)
