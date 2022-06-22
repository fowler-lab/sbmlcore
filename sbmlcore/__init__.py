#! /usr/bin/env python3

# from .AminoAcidProperties import *
from .AminoAcidProperties import AminoAcidProperty
from .AminoAcidProperties import AminoAcidVolumeChange
from .AminoAcidProperties import AminoAcidHydropathyChangeKyteDoolittle
from .AminoAcidProperties import AminoAcidHydropathyChangeWimleyWhite
from .AminoAcidProperties import AminoAcidMWChange
from .AminoAcidProperties import AminoAcidPiChange

from .ExternalCode import Stride
from .ExternalCode import FreeSASA
from .ExternalCode import SNAP2
from .StructuralDistances import StructuralDistances
from .TrajectoryDistances import TrajectoryDistances
from .TrajectoryDihedrals import TrajectoryDihedrals
from .DeepDDG import DeepDDG
from .FeaturesDataFrame import FeatureDataset

'''
Use of semantic versioning, MAJOR.MINOR.MAINTAINANCE where
MAJOR is not backwards compatible, but MINOR and MAINTAINANCE are
'''
__version__ = "0.0.1"
__author__ = 'Philip W Fowler and Charlotte I Lynch and Dylan Adlard'
