import pandas
import sbmlcore

def test_amino_acid_volume_change():

    a = {'mutation': ['A1D','E2K']}
    df = pandas.DataFrame.from_dict(a)

    a = sbmlcore.AminoAcidVolumeChange(df)
    df =  a.update()
    assert 'd_volume' in df.columns
