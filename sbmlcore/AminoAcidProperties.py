import pandas


class AminoAcidVolumeChange(object):

    def __init__(self, dataframe):

        assert isinstance(dataframe, pandas.DataFrame)
