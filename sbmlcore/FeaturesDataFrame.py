import pandas


class FeatureDataset(object):
    '''
    Dataset of protein mutations and associated structural- and chemical-features for simple machine learning.

    Args:
        dataframe (pandas.DataFrame) 
        protein (str, optional): protein name
        gene (str, optional): gene name
        species (str, optional): species
        reference (str)
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
        Overload the addition operator to allow sbml Features to be added.
        '''

        self.df = other._add_feature(self.df)

        return self

    def __repr__(self):
        '''
        Print a summary of this Features dataset to STDOUT.
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

        Args:
            sbmlcore.Feature or list of sbmlcore.Features
        '''

        if isinstance(other,list):
            for i in other:
                self.df = i._add_feature(self.df)
        else:
            self.df = other._add_feature(self.df)
