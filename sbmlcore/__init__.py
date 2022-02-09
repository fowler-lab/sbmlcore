#! /usr/bin/env python3

# from .AminoAcidProperties import *
from .AminoAcidProperties import AminoAcidProperty
from .AminoAcidProperties import AminoAcidVolumeChange
from .AminoAcidProperties import AminoAcidHydropathyChangeKyteDoolittle
from .AminoAcidProperties import AminoAcidHydropathyChangeWimleyWhite
from .AminoAcidProperties import AminoAcidMWChange
from .AminoAcidProperties import AminoAcidPiChange

from .ExternalCode import Stride
from .StructuralDistances import StructuralDistances

'''
Use of semantic versioning, MAJOR.MINOR.MAINTAINANCE where
MAJOR is not backwards compatible, but MINOR and MAINTAINANCE are
'''
__version__ = "0.0.1"
__author__ = 'Philip W Fowler and Charlotte I Lynch'
