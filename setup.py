#!/usr/bin/env python

from setuptools import setup

with open("README.md", "r") as f:
    README = f.read()

setup(
    name='sbmlcore',
    version='0.1.0',
    description='Infer if a non-resistant mutation is associated with a resistance mutation using Fishers Exact Test',
    author='Philip W Fowler and Charlotte I Lynch',
    url='https://github.com/fowler-lab/sbmlcore',
    long_description = README,
    install_requires=[
        'pandas',
        'pytest',
        'pytest-cov'
        ],
    packages = ['sbmlcore'],
    python_requires='>=3.8',
    zip_safe=False
    )
