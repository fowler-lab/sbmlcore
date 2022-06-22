import pandas

class FeatureDataset(object):
    '''
    Class that contains the dataset of mutations and associated features for simple machine learning.

    Parameters
    ----------
    existing dataframe

    Optional:
    protein (str)
    gene (str)
    species (str)
    reference (str)

    Returns
    -------
    dataframe with additional feature columns

    '''

    def __init__(self, df, *args, **kwargs):

        assert isinstance(df, pandas.DataFrame), 'must supply a pandas.DataFrame'

        assert 'mutation' in df.columns, 'passed pandas.DataFrame must contain a column called mutation'

        self.df = df

        allowed_kwargs = ['protein', 'gene', 'species', 'reference']

        seen = set()
        for key in kwargs.keys():
            if key in allowed_kwargs:
                setattr(self, key, kwargs[key])
                seen.add(key)

        for key in set(allowed_kwargs).difference(seen):
            # Default values to None if the kwarg has not been passed
            setattr(self, key, None)


    def __add__(self, other):
        ''''
        Overload the addition operator.
        '''
        self.df = other._add_feature(self.df)

        return self

    def __repr__(self):
        '''
        Overload the print function to write a summary of this Features dataset.

        Returns:
            str: String describing the Features dataset
        '''
        output = ''
        if self.species is not None:
            output += 'species:          ' + self.species + '\n'

        if self.gene is not None:
            output += 'gene name:        ' + self.gene + '\n'

        if self.protein is not None:
            output += 'protein name:     ' + self.protein + '\n'

        output += 'number of rows:   ' + str(len(self.df)) + '\n'
        output += '\n'

        if len(self.df) < 3:
            output += self.df.__repr__()
        else:
            output += self.df[:3].__repr__()
        return(output)


    def add_feature(self, other):
        '''
        Add the supplied sbmlcore Feature to this dataset.

        Arguments:
            sbmlcore.Feature or list of sbmlcore.Features
        '''
        if isinstance(other,list):
            for i in other:
                self.df = i._add_feature(self.df)
        else:
            self.df = other._add_feature(self.df)
