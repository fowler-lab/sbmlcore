import pandas

class AminoAcidProperty(object):

    def __init__(self, dataframe):

        assert isinstance(dataframe, pandas.DataFrame)

        self.dataframe = dataframe

        assert 'mutation' in self.dataframe.columns, 'passed dataframe must contain a column called mutations'

        def find_amino_acids(row):
            return(row.mutation[0], row.mutation[-1])

        self.dataframe[['ref_amino_acid','alt_amino_acid']] = self.dataframe.apply(find_amino_acids,axis=1)


class AminoAcidVolumeChange(AminoAcidProperty):

    def __init__(self):

        aa_volumes = {'A': 88.6, 'R': 173.4, 'N': 114.1, 'D': 111.1, 'C': 108.5,
                      'Q': 143.8, 'E': 138.4, 'G': 60.1, 'H': 153.2, 'I': 166.7,
                      'L': 166.7, 'K': 168.6, 'M': 162.9, 'F': 189.9, 'P': 112.7,
                      'S': 89.0, 'T': 116.1, 'W': 227.8, 'Y': 193.6, 'V': 140.0}

        rows = []
        for i in aa_volumes.keys():
            for j in aa_volumes.keys():
                rows.append([i, j, aa_volumes[i] - aa_volumes[j]])

        self.lookup = pandas.DataFrame(rows,columns=['ref_amino_acid', 'alt_amino_acid', 'd_volume'])

        self.lookup.set_index(['ref_amino_acid', 'alt_amino_acid'],inplace=True)


    def __add__(self, other):

        assert isinstance(other, pandas.DataFrame)

        assert 'mutation' in other.columns, 'passed dataframe must contain a column called mutations'

        def find_amino_acids(row):
            return(row.mutation[0], row.mutation[-1])

        other[['ref_amino_acid','alt_amino_acid']] = other.apply(find_amino_acids,axis=1)

        other.set_index(['ref_amino_acid', 'alt_amino_acid'], inplace=True)

        other = other.join(self.lookup)

        other.reset_index(inplace=True)

        other.drop(columns = ['ref_amino_acid', 'alt_amino_acid'], inplace=True)

        return(other)

    def __radd__(self, other):

        return self.__add__(other)
