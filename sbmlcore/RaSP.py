import pathlib
import copy

import pandas

class RaSP(object):
    """
    Prediction of protein stability via the RaSP Google Colab project.

    Args:
        DeepDDGFile (file): file written out by the DeepDDG webserver
        offsets (dict): dictionary of form {segid (str): value (int)} where value is 
                the numerical offset between the genetic sequence and the PDB
    """

    def __init__(self, RaSPFile, offsets=None):

        assert pathlib.Path(RaSPFile).is_file(), "specified file does not exist!"

        self.rasp_file = RaSPFile

        self.results = pandas.read_csv(self.rasp_file)

        self.results.drop(columns=['pdbid', 'variant', 'wt_idx', 'mt_idx', 'wt'], inplace=True)

        self.results.rename(columns={'chainid': 'segid',\
                                     'pos': 'resid',\
                                     'wt_AA': 'ref_amino_acid',\
                                     'mt_AA': 'alt_amino_acid',\
                                     'score_ml_fermi': 'rasp_score_ml_fermi',\
                                     'score_ml': 'rasp_score_ml',\
                                     'wt_nlf': 'rasp_wt_nlf',\
                                     'mt_nlf': 'rasp_mt_nlf'}, inplace=True)
       
        # apply any offsets to the residue numbering
        # as specified in the supplied offsets dict e.g. {'A': 3, 'B': -4}
        # Chain is segid i.e. A, B, C etc.

        def update_resid(row, offsets):
            offset = offsets[row['segid']]
            return row['resid'] + offset

        if offsets is not None:

            assert isinstance(offsets, dict), "Offsets should be specified as a dictionary e.g. offsets = {'A': 3, 'B': -4}"

            assert set(offsets.keys()) == set(self.results.segid.unique()), 'specified segids in offsets do not match what is in the provided RaSP file!'

            for chain in offsets:
                assert isinstance(offsets[chain], int), "Offsets for each segid must be an integer!"

            self.results['resid'] = self.results.apply(update_resid, args=(offsets,), axis=1)

    def _add_feature(self, other):

        assert isinstance(other, pandas.DataFrame)

        assert 'mutation' in other.columns, 'passed dataframe must contain a column called mutation'

        assert 'segid' in other.columns, 'passed dataframe must contain a column called segid containing chain information e.g. A'

        assert 'rasp_score_ml' not in other.columns, 'passed dataframe already contains a RaSP column -> have you already added it?'

        def split_mutation(row):
            return pandas.Series([row.mutation[0], int(row.mutation[1:-1]), row.mutation[-1]])

        other[['ref_amino_acid', 'resid', 'alt_amino_acid']] = other.apply(split_mutation, axis=1)
        other.set_index(['segid', 'ref_amino_acid', 'resid', 'alt_amino_acid'], inplace=True)
        
        self.results.set_index(['segid', 'ref_amino_acid', 'resid', 'alt_amino_acid'], inplace=True)

        other = other.join(self.results, how='left')

        other.reset_index(inplace=True)
        self.results.reset_index(inplace=True)

        other.drop(columns = ['ref_amino_acid', 'resid', 'alt_amino_acid'], inplace=True)

        return(other)
