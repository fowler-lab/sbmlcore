import pathlib
import copy

import pandas

class DeepDDG(object):
    """
    Prediction of protein stability.

    Parameters
    ----------
    .ddg file
    offsets - as dictionary of form {segid (str): value (int)}

    Returns
    -------
    Dataframe with protein stability as an extra column

    """

    def __init__(self, DeepDDGFile, offsets=None):

        assert pathlib.Path(DeepDDGFile).is_file(), "specified file does not exist!"

        self.deep_ddg_file = DeepDDGFile

        # unfortunately deepDDG uses a non-standard tab/space delimited file format
        # this breaks pandas.read_csv so it is easier to do it the long way....
        results = []

        # open the file for reading
        with open(DeepDDGFile, 'r') as INPUT:

            # read the first line and check it contains what we expect
            header = INPUT.readline()
            assert "#chain WT ResID Mut ddG" in header, 'file does not contain expected columns -> are you sure this is a DeepDDG output file?'

            # iterate through the remainder of the file
            for i in INPUT:

                cols = i.rstrip().split(' ')

                # the first few columns always line up...
                line = [cols[0], cols[1], int(cols[2]), cols[3]]

                # ..but negative values of ddG can be found in the 7th column
                if len(cols) == 6:
                    line.append(float(cols[5]))
                # ..and positive the 8th
                elif len(cols) == 7:
                    line.append(float(cols[6]))
                else:
                    raise IOError("DeepDDG file not formatted as expected!")

                results.append(line)

        # finally build and store the dataframe
        self.results =pandas.DataFrame(results, columns=['segid', 'ref_amino_acid', 'resid', 'alt_amino_acid', 'deep_ddG'])

        # apply any offsets to the residue numbering
        # as specified in the supplied offsets dict e.g. {'A': 3, 'B': -4}
        # Chain is segid i.e. A, B, C etc.

        def update_resid(row, offsets):
            offset = offsets[row['segid']]
            return row['resid'] + offset

        if offsets is not None:

            assert isinstance(offsets, dict), "Offsets should be specified as a dictionary e.g. offsets = {'A': 3, 'B': -4}"

            assert set(offsets.keys()) == set(self.results.segid.unique()), 'specified segids in offsets do not match what is in the provided DeepDDG file!'

            for chain in offsets:
                assert isinstance(offsets[chain], int), "Offsets for each segid must be an integer!"

            self.results['resid'] = self.results.apply(update_resid, args=(offsets,), axis=1)

    def _add_feature(self, other):

        assert isinstance(other, pandas.DataFrame)

        assert 'mutation' in other.columns, 'passed dataframe must contain a column called mutation'

        assert 'segid' in other.columns, 'passed dataframe must contain a column called segid containing chain information e.g. A'

        assert 'deep_ddG' not in other.columns, 'passed dataframe already contains a deep_ddG column -> have you already added it?'

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
