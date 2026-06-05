import importlib
import sys

import pandas


def test_top_level_import_exposes_lightweight_features():
    sys.modules.pop("sbmlcore", None)

    sbmlcore = importlib.import_module("sbmlcore")

    feature = sbmlcore.AminoAcidMWChange()
    df = pandas.DataFrame({"mutation": ["A1D"]})
    result = feature._add_feature(df)

    assert "d_MW" in result.columns
    assert result.loc[0, "d_MW"] == 44.0


def test_residue_depth_accepts_offsets_without_segids(monkeypatch):
    module = importlib.import_module("sbmlcore.ResidueDepth")

    class FakePath:
        def is_file(self):
            return True

    class FakeResidue:
        def __init__(self, resid):
            self.id = (" ", resid, " ")

    class FakeChain:
        def __init__(self, segid, resids):
            self._segid = segid
            self._residues = [FakeResidue(resid) for resid in resids]

        def get_id(self):
            return self._segid

        def get_residues(self):
            return list(self._residues)

        def __getitem__(self, resid):
            return FakeResidue(resid)

    class FakeModel:
        def __init__(self):
            self._chains = {"A": FakeChain("A", [1, 2])}

        def get_chains(self):
            return list(self._chains.values())

        def __getitem__(self, segid):
            return self._chains[segid]

    class FakeStructure:
        def __getitem__(self, index):
            assert index == 0
            return FakeModel()

    class FakeParser:
        def get_structure(self, *_args, **_kwargs):
            return FakeStructure()

    monkeypatch.setattr(module.pathlib, "Path", lambda *_args, **_kwargs: FakePath())
    monkeypatch.setattr(module, "PDBParser", lambda: FakeParser())
    monkeypatch.setattr(module, "get_surface", lambda _model: object())
    monkeypatch.setattr(
        module, "residue_depth", lambda residue, _surface: residue.id[1] / 10
    )

    result = module.ResidueDepth("fake.pdb", offsets={"A": 10}).results

    assert result.to_dict(orient="records") == [
        {"segid": "A", "resid": 11, "depth": 0.1},
        {"segid": "A", "resid": 12, "depth": 0.2},
    ]
