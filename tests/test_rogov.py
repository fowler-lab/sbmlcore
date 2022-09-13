import pandas
import numpy
import pytest
import sbmlcore

def test_amino_acid_rogov_change_value():

    a = {'mutation': ['A1D', 'E2K']}
    df = pandas.DataFrame(a)

    a = sbmlcore.AminoAcidRogovChange()
    df =  a._add_feature(df)
    assert 'd_rogov' in df.columns

    assert df.d_rogov.values == pytest.approx(numpy.array([0.029, 0.224]))

