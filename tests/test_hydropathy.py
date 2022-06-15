# Tests for the '.AminoAcidHydropathyChange' defs for different hydropathy scales
import pandas
import numpy
import pytest
import sbmlcore

def test_amino_acid_hydropathy_KD_change_value():

    a = {'mutation': ['A1D', 'E2K']}
    df = pandas.DataFrame(a)

    a = sbmlcore.AminoAcidHydropathyChangeKyteDoolittle()
    df =  a._add_feature(df)
    assert 'd_hydropathy_KD' in df.columns

    assert df.d_hydropathy_KD.values == pytest.approx(numpy.array([-5.3, -0.4]))

    # these should all fail!
    with pytest.raises(AssertionError):
        assert df.d_hydropathy_KD.values == pytest.approx(numpy.array([5.3, 0.4]))

    with pytest.raises(AssertionError):
        assert df.d_hydropathy_KD.values == pytest.approx(numpy.array([5.3, -0.4]))

    with pytest.raises(AssertionError):
        assert df.d_hydropathy_KD.values == pytest.approx(numpy.array([-5.3, '-0.4']))

def test_amino_acid_hydropathy_WW_change_value():

    a = {'mutation': ['A1D', 'E2K']}
    df = pandas.DataFrame(a)

    a = sbmlcore.AminoAcidHydropathyChangeWimleyWhite()
    df =  a._add_feature(df)
    assert 'd_hydropathy_WW' in df.columns

    assert df.d_hydropathy_WW.values == pytest.approx(numpy.array([2.08, 0.20]))

    # these should all fail!
    with pytest.raises(AssertionError):
        assert df.d_hydropathy_WW.values == pytest.approx(numpy.array([-2.08, -0.20]))

    with pytest.raises(AssertionError):
        assert df.d_hydropathy_WW.values == pytest.approx(numpy.array([2.08, -0.20]))

    with pytest.raises(AssertionError):
        assert df.d_hydropathy_WW.values == pytest.approx(numpy.array([2.08, '0.20']))
