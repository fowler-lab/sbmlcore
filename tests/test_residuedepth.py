import pandas
import numpy
import pytest
import sbmlcore

def test_residue_depth():

    a = sbmlcore.ResidueDepth('tests/3pl1-mod-chains.pdb', segids=['A'])

    assert float(a.results[a.results.resid==1].depth) == pytest.approx(2.10304859)

    