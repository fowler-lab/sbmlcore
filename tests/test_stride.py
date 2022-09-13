import pytest

import pandas

import sbmlcore


def test_stride():

    a = sbmlcore.Stride('tests/3pl1.pdb')

    assert len(a.results) == 185

    # check that all the -180 < phi <= 180 (and also for psi)
    assert all((a.results.phi > -180) & (a.results.phi <= 180))

    assert (a.results.n_hbond_acceptors>=0).sum() == 185
    assert (a.results.n_hbond_donors>=0).sum() == 185

    # check that only these 1-letter codes are used for secondary structure
    assert set(a.results.secondary_structure).issubset({'H', 'I', 'G', 'E', 'B', 'T', 'C'})

    # these should all fail!
    with pytest.raises(AssertionError):
        a = sbmlcore.Stride('tests/3pl2.pdb')

    b = {'segid': ['A', 'A', 'A'], 'mutation': ['M1D','R2K', 'A3V']}
    df = pandas.DataFrame.from_dict(b)

    # test out individual features
    df = a._add_feature(df, 'psi')
    assert 'psi' in df.columns

    df = a._add_feature(df, 'phi')
    assert 'phi' in df.columns

    b = {'segid': ['A', 'A', 'A'], 'mutation': ['M1D','R2K', 'A3V']}
    df = pandas.DataFrame.from_dict(b)

    # test adding all
    df = a._add_feature(df)
    assert 'psi' in df.columns
    assert 'phi' in df.columns
