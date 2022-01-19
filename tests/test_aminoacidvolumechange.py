import pandas
import numpy
import pytest
import sbmlcore

def test_amino_acid_volume_change_df_format():

    a = {'mutation': ['A1D','E2K']}
    df = pandas.DataFrame.from_dict(a)

    a = sbmlcore.AminoAcidVolumeChange()
    df2 =  a + df
    assert 'd_volume' in df2.columns

def test_amino_acid_volume_change_value():

    a = {'mutation': ['A1D', 'E2K']}
    df = pandas.DataFrame.from_dict(a)

    a = sbmlcore.AminoAcidVolumeChange()
    df2 = a + df
    assert 'd_volume' in df2.columns

    assert df2.d_volume.values == pytest.approx(numpy.array([-22.5, -30.2]))

    # this should fail!
    with pytest.raises(AssertionError):
        assert df2.d_volume.values == pytest.approx(numpy.array([22.5, 30.2]))
