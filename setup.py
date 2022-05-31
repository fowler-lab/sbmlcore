#!/usr/bin/env python

from setuptools import setup

with open("README.md", "r") as f:
    README = f.read()

setup(
    name='sbmlcore',
    version='0.1.0',
    description='Constructs core features table for the application to machine learning models',
    author='Philip W Fowler and Charlotte I Lynch and Dylan Adlard',
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
