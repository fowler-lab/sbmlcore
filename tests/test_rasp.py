import pytest

import pandas

import sbmlcore


def test_rasp():

    a = sbmlcore.RaSP('tests/cavity_pred_3PL1_A.csv')

    assert len(a.results) == 185*20

    assert all((a.results.rasp_score_ml_fermi >= 0) & (a.results.rasp_score_ml_fermi <= 1))

    assert all(a.results.rasp_score_ml >= -2)    

    # these should all fail!
    with pytest.raises(AssertionError):
        a = sbmlcore.RaSP('tests/3pl2.pdb')

    b = {'segid': ['A', 'A', 'A'], 'mutation': ['M1D','R2K', 'A3V']}
    df = pandas.DataFrame.from_dict(b)

    # test out individual features
    df = a._add_feature(df)
    assert 'rasp_score_ml' in df.columns
    assert 'rasp_score_ml_fermi' in df.columns

