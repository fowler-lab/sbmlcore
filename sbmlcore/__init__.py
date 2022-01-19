#! /usr/bin/env python3
class sbmlcore:
    """Overall module framework for mutation prediction."""

    name: str = "Structural biology machine learning core framework"

    def __init__(self):
        self.sbmlcore = sbmlcore

from .AminoAcidProperties import AminoAcidProperty
from .AminoAcidProperties import AminoAcidVolumeChange

'''
Use of semantic versioning, MAJOR.MINOR.MAINTAINANCE where
MAJOR is not backwards compatible, but MINOR and MAINTAINANCE are
'''
__version__ = "0.0.1"
__author__ = 'Philip W Fowler'
