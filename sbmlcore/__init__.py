#! /usr/bin/env python3

from importlib import import_module
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

__author__ = "Philip W Fowler and Charlotte I Lynch and Dylan Adlard"

_EXPORTS = {
    "AminoAcidPropertyChange": "sbmlcore.AminoAcidProperties",
    "AminoAcidVolumeChange": "sbmlcore.AminoAcidProperties",
    "AminoAcidHydropathyChangeKyteDoolittle": "sbmlcore.AminoAcidProperties",
    "AminoAcidHydropathyChangeWimleyWhite": "sbmlcore.AminoAcidProperties",
    "AminoAcidMWChange": "sbmlcore.AminoAcidProperties",
    "AminoAcidPiChange": "sbmlcore.AminoAcidProperties",
    "AminoAcidRogovChange": "sbmlcore.AminoAcidProperties",
    "AminoAcidVolume": "sbmlcore.AminoAcidProperties",
    "AminoAcidHydropathyKyteDoolittle": "sbmlcore.AminoAcidProperties",
    "AminoAcidHydropathyWimleyWhite": "sbmlcore.AminoAcidProperties",
    "AminoAcidMW": "sbmlcore.AminoAcidProperties",
    "AminoAcidPi": "sbmlcore.AminoAcidProperties",
    "SideChainRings": "sbmlcore.AminoAcidProperties",
    "HBondDonors": "sbmlcore.AminoAcidProperties",
    "HBondAcceptors": "sbmlcore.AminoAcidProperties",
    "amino_acid_3to1letter": "sbmlcore.Misc",
    "amino_acid_1to3letter": "sbmlcore.Misc",
    "Stride": "sbmlcore.ExternalCode",
    "FreeSASA": "sbmlcore.ExternalCode",
    "SNAP2": "sbmlcore.ExternalCode",
    "TempFactors": "sbmlcore.TempFactors",
    "StructuralDistances": "sbmlcore.StructuralDistances",
    "TrajectoryDistances": "sbmlcore.TrajectoryDistances",
    "TrajectoryDihedrals": "sbmlcore.TrajectoryDihedrals",
    "DeepDDG": "sbmlcore.DeepDDG",
    "RaSP": "sbmlcore.RaSP",
    "ResidueDepth": "sbmlcore.ResidueDepth",
    "FeatureDataset": "sbmlcore.FeaturesDataFrame",
}

__all__ = sorted(_EXPORTS)


def _read_local_version():
    version_file = Path(__file__).resolve().parent.parent / "VERSION"
    if not version_file.exists():
        raise FileNotFoundError(version_file)

    raw_version = version_file.read_text(encoding="utf-8").strip()
    return raw_version.lstrip("v")


try:
    __version__ = _read_local_version()
except FileNotFoundError:
    try:
        __version__ = version("sbmlcore")
    except PackageNotFoundError:
        __version__ = "0.0.0"


def __getattr__(name):
    if name not in _EXPORTS:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

    module = import_module(_EXPORTS[name])
    value = getattr(module, name)
    globals()[name] = value
    return value


def __dir__():
    return sorted(list(globals().keys()) + __all__)
