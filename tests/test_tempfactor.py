import pytest

import pandas

import sbmlcore


def test_tempfactors():

    a = sbmlcore.TempFactors('tests/3pl1.pdb')

    assert len(a.results) == 185

    # check that all the -180 < phi <= 180 (and also for psi)
    assert all(a.results.temp_factor >= 0)

    assert float(a.results[a.results.resid==1].temp_factor) == pytest.approx(59.72)

    # these should all fail!
    with pytest.raises(AssertionError):
        a = sbmlcore.TempFactors('tests/3pl2.pdb')

    b = {'segid': ['A', 'A', 'A'], 'mutation': ['M1D','R2K', 'A3V']}
    df = pandas.DataFrame.from_dict(b)

    # test adding all
    df = a._add_feature(df)
    assert 'temp_factor' in df.columns
