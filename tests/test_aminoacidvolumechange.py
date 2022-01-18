import pandas
from sbmlcore.AminoAcidProperties import AminoAcidProperty

def test_amino_acid_volume_change_df_format():

    a = {'mutation': ['A1D','E2K']}
    df = pandas.DataFrame.from_dict(a)

    a = sbmlcore.AminoAcidVolumeChange()
    df3 =  a + df
    assert 'd_volume' in df3.columns

def test_amino_acid_volume_change_value():

    a = {'mutation': ['A1D']}
    df = pandas.Dataframe.from_dict(a)

    a = sbmlcore.AminoAcidVolumeChange()
    df2 = a + df
    assert 'd_volume' in df2.columns == 22.5
    #Currently D-A, but AminoAcidProperties uses A-D therefore this test should fail. 
