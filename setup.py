# setup.py

from setuptools import setup, find_packages

setup(
    name='hypermatrix',
    version='0.1',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'hypermatrix = hypermatrix.main:main',
        ],
    },
    install_requires=[
        # List any dependencies

    ],
)

