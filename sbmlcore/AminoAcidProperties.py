import pkg_resources

import pandas
#N.B. All property scales have been checked! CIL

def _make_table(data_dict):
    rows = []
    for i in data_dict.keys():
        for j in data_dict.keys():
            rows.append([i, j, data_dict[j] - data_dict[i]])
    return pandas.DataFrame(rows, columns=['ref_amino_acid', 'alt_amino_acid', 'value'])


class AminoAcidProperty(object):
    """Base class for storing amino acid properties.
    
    Notes:
        uses a pandas.DataFrame
    """

    def __init__(self):

        pass

    def _add_feature(self, other):

        assert isinstance(other, pandas.DataFrame)

        assert self.lookup.columns[0] not in other.columns, 'trying to add a column that is already in the dataset!'

        assert 'mutation' in other.columns, 'passed dataframe must contain a column called mutation'

        def find_amino_acids(row):
            return(pandas.Series([row.mutation[0], row.mutation[-1]]))

        other[['ref_amino_acid','alt_amino_acid']] = other.apply(find_amino_acids,axis=1)

        other.set_index(['ref_amino_acid', 'alt_amino_acid'], inplace=True)

        other = other.join(self.lookup, how='left')

        other.reset_index(inplace=True)

        other.drop(columns = ['ref_amino_acid', 'alt_amino_acid'], inplace=True)

        return(other)


class AminoAcidRogovChange(AminoAcidProperty):

    def __init__(self):

        filename = pkg_resources.resource_filename("sbmlcore", 'data/rogov.csv')
        self.lookup = pandas.read_csv(filename)

        def split_row(row):
            if isinstance(row.MUTATION, str) and len(row.MUTATION)==2:
                return(pandas.Series([row.MUTATION[0], row.MUTATION[-1]]))
            else:
                return(pandas.Series([None,None]))


        self.lookup[['ref_amino_acid', 'alt_amino_acid']] = self.lookup.apply(split_row, axis=1)
        self.lookup.rename(columns={'SCORE': 'd_rogov'}, inplace=True)
        self.lookup.drop(columns=['MUTATION'], inplace=True)
        self.lookup.set_index(['ref_amino_acid', 'alt_amino_acid'], inplace=True)


class AminoAcidVolumeChange(AminoAcidProperty):
    """
    Change in sidechain volume (ang^3) of an amino acid due to a mutation.
    """

    def __init__(self):

        aa_volume = {'A': 88.6, 'R': 173.4, 'N': 114.1, 'D': 111.1, 'C': 108.5,
                      'Q': 143.8, 'E': 138.4, 'G': 60.1, 'H': 153.2, 'I': 166.7,
                      'L': 166.7, 'K': 168.6, 'M': 162.9, 'F': 189.9, 'P': 112.7,
                      'S': 89.0, 'T': 116.1, 'W': 227.8, 'Y': 193.6, 'V': 140.0}

        self.lookup = _make_table(aa_volume)
        self.lookup.rename(columns = {'value': 'd_volume'}, inplace=True)
        self.lookup.set_index(['ref_amino_acid', 'alt_amino_acid'],inplace=True)


class AminoAcidHydropathyChangeKyteDoolittle(AminoAcidProperty):
    """
    Change in hydropathy due to a mutation.

    Notes:
        Kyte and Doolittle (1982) DOI 10.1016/0022-2836(82)90515-0
        Wimley-White is also available.
    """

    def __init__(self):

        aa_hydropathy_KD = {'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5,
                            'C': 2.5, 'E': -3.5, 'Q': -3.5, 'G': -0.4,
                            'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9,
                            'M': 1.9, 'F': 2.8, 'P': -1.6, 'S': -0.8,
                            'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2}

        self.lookup = _make_table(aa_hydropathy_KD)
        self.lookup.rename(columns = {'value': 'd_hydropathy_KD'}, inplace=True)
        self.lookup.set_index(['ref_amino_acid', 'alt_amino_acid'],inplace=True)


class AminoAcidHydropathyChangeWimleyWhite(AminoAcidProperty):
    """
    Change in hydropathy of an amino acid site due to a mutation.

    Notes:
        Wimley White (1996) DOI 10.1038/nsb1096-842, octanol-interface scale from https://blanco.biomol.uci.edu/hydrophobicity_scales.html, Asp- Glu- His+ used.
        Kyte and Doolittle (1982) is also available.
    """

    def __init__(self):

        aa_hydropathy_WW = {'A': 0.33, 'R': 1.00, 'N': 0.43, 'D': 2.41,
                            'C': 0.22, 'Q': 0.19, 'E': 1.61, 'G': 1.14,
                            'H': 1.37, 'I': -0.81, 'L': -0.69, 'K': 1.81,
                            'M': -0.44, 'F': -0.58, 'P': -0.31, 'S': 0.33,
                            'T': 0.11, 'W': -0.24, 'Y': 0.23, 'V': -0.53}

        self.lookup = _make_table(aa_hydropathy_WW)
        self.lookup.rename(columns = {'value': 'd_hydropathy_WW'}, inplace=True)
        self.lookup.set_index(['ref_amino_acid', 'alt_amino_acid'],inplace=True)


class AminoAcidMWChange(AminoAcidProperty):
    """
    Change in molecular weight (g/mol) of an amino acid site due to a mutation.

    Notes:
        From: https://www.thermofisher.com/uk/en/home/references/ambion-tech-support/rna-tools-and-calculators/proteins-and-amino-acids.html
    """

    def __init__(self):

        aa_MW = {'A': 89.1, 'R': 174.2, 'N': 132.1, 'D': 133.1, 'C': 121.2,
                 'E': 147.1, 'Q': 146.2, 'G':75.1, 'H': 155.2, 'I': 131.2,
                 'L': 131.2, 'K': 146.2, 'M': 149.2, 'F': 165.2, 'P': 115.1,
                 'S': 105.1, 'T': 119.1, 'W': 204.2, 'Y': 181.2, 'V': 117.1}

        self.lookup = _make_table(aa_MW)
        self.lookup.rename(columns = {'value': 'd_MW'}, inplace=True)
        self.lookup.set_index(['ref_amino_acid', 'alt_amino_acid'],inplace=True)


class AminoAcidPiChange(AminoAcidProperty):
    """
    Change in isoelectric point of an amino acid site due to a mutation.

    Notes:
        From https://www.sigmaaldrich.com/life-science/metabolomics/learning-center/amino-acid-reference-chart.html
    """

    def __init__(self):

        aa_Pi = {'A': 6.00, 'R': 10.76, 'N': 5.41, 'D': 2.77, 'C': 5.07, 'E': 3.22,
             'Q': 5.65, 'G': 5.97, 'H': 7.59, 'I': 6.02, 'L': 5.98, 'K': 9.74,
             'M': 5.74, 'F': 5.48, 'P': 6.30, 'S': 5.68, 'T': 5.60, 'W': 5.89,
             'Y': 5.66, 'V': 5.96}

        self.lookup = _make_table(aa_Pi)
        self.lookup.rename(columns = {'value': 'd_Pi'}, inplace=True)
        self.lookup.set_index(['ref_amino_acid', 'alt_amino_acid'],inplace=True)
