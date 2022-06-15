# Tests for the '.AminoAcidPiChange' defs
import pandas
import numpy
import pytest
import sbmlcore

def test_amino_acid_Pi_change_value():

    a = {'mutation': ['A1D', 'E2K']}
    df = pandas.DataFrame(a)

    a = sbmlcore.AminoAcidPiChange()
    df =  a._add_feature(df)
    assert 'd_Pi' in df.columns

    assert df.d_Pi.values == pytest.approx(numpy.array([-3.23, 6.52]))

    # these should all fail!
    with pytest.raises(AssertionError):
        assert df.d_Pi.values == pytest.approx(numpy.array([3.23, -6.52]))

    with pytest.raises(AssertionError):
        assert df.d_Pi.values == pytest.approx(numpy.array([-3.23, -6.52]))

    with pytest.raises(AssertionError):
        assert df.d_Pi.values == pytest.approx(numpy.array([-3.23, '6.52']))
