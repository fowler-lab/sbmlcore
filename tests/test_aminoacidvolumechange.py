import smblcore

def test_amino_acid_volume_change():

    a = {'MUTATIONS': ['A1D','E2K']}
    df = pandas.DataFrame.from_dict(a)

    df = smblcore.AminoAcidVolumeChange(df)

    assert ['d_volume'] in df.columns
    
