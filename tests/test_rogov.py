import pandas
import numpy
import pytest
import sbmlcore

def test_amino_acid_roskov_change_value():

    a = {'mutation': ['A1D', 'E2K']}
    df = pandas.DataFrame(a)

    a = sbmlcore.AminoAcidRoskovChange()
    df =  a._add_feature(df)
    assert 'd_noskov' in df.columns

    assert df.d_noskov.values == pytest.approx(numpy.array([0.029, 0.224]))

