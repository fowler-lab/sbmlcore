# Tests for the '.AminoAcidMWChange' defs
import pandas
import numpy
import pytest
import sbmlcore

def test_amino_acid_MW_change_value():

    a = {'mutation': ['A1D', 'E2K']}
    df = pandas.DataFrame(a)

    a = sbmlcore.AminoAcidMWChange()
    df =  a._add_feature(df)
    assert 'd_MW' in df.columns

    assert df.d_MW.values == pytest.approx(numpy.array([44.0, -0.9]))

    # these should all fail!
    with pytest.raises(AssertionError):
        assert df.d_MW.values == pytest.approx(numpy.array([-44.0, 0.9]))

    with pytest.raises(AssertionError):
        assert df.d_MW.values == pytest.approx(numpy.array([-44.0, -0.9]))

    with pytest.raises(AssertionError):
        assert df.d_MW.values == pytest.approx(numpy.array([44.0, '-0.9']))
