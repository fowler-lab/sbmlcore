# Tests for .AminoAcidVolumeChange
import pandas
import numpy
import pytest
import sbmlcore

def test_amino_acid_volume_change_value():

    a = {'mutation': ['A1D', 'E2K']}
    df = pandas.DataFrame.from_dict(a)

    a = sbmlcore.AminoAcidVolumeChange()
    df =  a._add_feature(df)
    assert 'd_volume' in df.columns

    assert df.d_volume.values == pytest.approx(numpy.array([22.5, 30.2]))

    # these should all fail!
    with pytest.raises(AssertionError):
        assert df.d_volume.values == pytest.approx(numpy.array([-22.5, -30.2]))

    with pytest.raises(AssertionError):
        assert df.d_volume.values == pytest.approx(numpy.array([22.5, -30.2]))

    with pytest.raises(AssertionError):
        assert df.d_volume.values == pytest.approx(numpy.array([22.5, '30.2']))
