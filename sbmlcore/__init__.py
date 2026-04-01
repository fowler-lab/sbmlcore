#! /usr/bin/env python3

from importlib import import_module
from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as package_version
from pathlib import Path


def _resolve_version():
    version_file = Path(__file__).resolve().parent.parent / "VERSION"
    if version_file.is_file():
        return version_file.read_text(encoding="utf-8").strip()

    try:
        return package_version("sbmlcore")
    except PackageNotFoundError:
        return "0+unknown"

_EXPORTS = {
    "AminoAcidPropertyChange": (".AminoAcidProperties", "AminoAcidPropertyChange"),
    "AminoAcidVolumeChange": (".AminoAcidProperties", "AminoAcidVolumeChange"),
    "AminoAcidHydropathyChangeKyteDoolittle": (
        ".AminoAcidProperties",
        "AminoAcidHydropathyChangeKyteDoolittle",
    ),
    "AminoAcidHydropathyChangeWimleyWhite": (
        ".AminoAcidProperties",
        "AminoAcidHydropathyChangeWimleyWhite",
    ),
    "AminoAcidMWChange": (".AminoAcidProperties", "AminoAcidMWChange"),
    "AminoAcidPiChange": (".AminoAcidProperties", "AminoAcidPiChange"),
    "AminoAcidRogovChange": (".AminoAcidProperties", "AminoAcidRogovChange"),
    "AminoAcidVolume": (".AminoAcidProperties", "AminoAcidVolume"),
    "AminoAcidHydropathyKyteDoolittle": (
        ".AminoAcidProperties",
        "AminoAcidHydropathyKyteDoolittle",
    ),
    "AminoAcidHydropathyWimleyWhite": (
        ".AminoAcidProperties",
        "AminoAcidHydropathyWimleyWhite",
    ),
    "AminoAcidMW": (".AminoAcidProperties", "AminoAcidMW"),
    "AminoAcidPi": (".AminoAcidProperties", "AminoAcidPi"),
    "SideChainRings": (".AminoAcidProperties", "SideChainRings"),
    "HBondDonors": (".AminoAcidProperties", "HBondDonors"),
    "HBondAcceptors": (".AminoAcidProperties", "HBondAcceptors"),
    "amino_acid_3to1letter": (".Misc", "amino_acid_3to1letter"),
    "amino_acid_1to3letter": (".Misc", "amino_acid_1to3letter"),
    "Stride": (".ExternalCode", "Stride"),
    "FreeSASA": (".ExternalCode", "FreeSASA"),
    "SNAP2": (".ExternalCode", "SNAP2"),
    "TempFactors": (".TempFactors", "TempFactors"),
    "StructuralDistances": (".StructuralDistances", "StructuralDistances"),
    "TrajectoryDistances": (".TrajectoryDistances", "TrajectoryDistances"),
    "TrajectoryDihedrals": (".TrajectoryDihedrals", "TrajectoryDihedrals"),
    "DeepDDG": (".DeepDDG", "DeepDDG"),
    "RaSP": (".RaSP", "RaSP"),
    "ResidueDepth": (".ResidueDepth", "ResidueDepth"),
    "FeatureDataset": (".FeaturesDataFrame", "FeatureDataset"),
}


def __getattr__(name):
    if name not in _EXPORTS:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

    module_name, attribute_name = _EXPORTS[name]
    value = getattr(import_module(module_name, __name__), attribute_name)
    globals()[name] = value
    return value


__all__ = sorted(_EXPORTS)

'''
Use of semantic versioning, MAJOR.MINOR.MAINTAINANCE where
MAJOR is not backwards compatible, but MINOR and MAINTAINANCE are
'''
__version__ = _resolve_version()
__author__ = 'Philip W Fowler and Charlotte I Lynch and Dylan Adlard'
